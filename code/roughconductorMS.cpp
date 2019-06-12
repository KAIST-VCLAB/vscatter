/****************************************************************************
- Codename: Multiple Scattering in V-groove (SIGGRAPH Asia 2018)

- Writers:   Joo Ho Lee(jhlee@vclab.kaist.ac.kr), Min H. Kim (minhkim@vclab.kaist.ac.kr)

- Institute: KAIST Visual Computing Laboratory

- Bibtex:
@Article{V-Scattering:SIGA:2018,
  author  = {Joo Ho Lee and Adrian Jarabo and Daniel S. Jeon and Diego Gutierrez and Min H. Kim},
  title   = {Practical Multiple Scattering for Rough Surfaces},
  journal = {ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2018)},
  year    = {2018},
  volume  = {36},
  number  = {6},
  pages   = {275:1--12},
  doi     = "10.1145/3272127.3275016",
  url     = "http://dx.doi.org/10.1145/3272127.3275016",
} 

- Joo Ho Lee and Min H. Kim have developed this software and related documentation
  (the "Software"); confidential use in source form of the Software,
  without modification, is permitted provided that the following
  conditions are met:
  1. Neither the name of the copyright holder nor the names of any
  contributors may be used to endorse or promote products derived from
  the Software without specific prior written permission.
  2. The use of the software is for Non-Commercial Purposes only. As
  used in this Agreement, "Non-Commercial Purpose" means for the
  purpose of education or research in a non-commercial organisation
  only. "Non-Commercial Purpose" excludes, without limitation, any use
  of the Software for, as part of, or in any way in connection with a
  product (including software) or service which is sold, offered for
  sale, licensed, leased, published, loaned or rented. If you require
  a license for a use excluded by this agreement,
  please email [minhkim@kaist.ac.kr].
  
- License:  GNU General Public License Usage
  Alternatively, this file may be used under the terms of the GNU General
  Public License version 3.0 as published by the Free Software Foundation
  and appearing in the file LICENSE.GPL included in the packaging of this
  file. Please review the following information to ensure the GNU General
  Public License version 3.0 requirements will be met:
  http://www.gnu.org/copyleft/gpl.html.

- Warranty: KAIST-VCLAB MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE 
  SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT 
  LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
  PARTICULAR PURPOSE, OR NON-INFRINGEMENT. KAIST-VCLAB SHALL NOT BE LIABLE FOR ANY 
  DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
  THIS SOFTWARE OR ITS DERIVATIVES
*****************************************************************************/

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

#include "vscatter.h"

static inline Spectrum computeMultipleFresnel(float cos_phi_sp, float sin_phi_sp, float theta_s, vec3 wi, Spectrum eta, Spectrum k, int bounce_n, bool isLeft) {

	theta_s = abs(theta_s);

	if (isLeft) {

		//required variable
		//cos_phi_sp
		//sin_phi_sp
		//theta_s
		//s0
		Spectrum F(1.f);

		float theta_v = M_PI - 2 * abs(theta_s);

		float cos_theta_v = cos(theta_v);
		float sin_theta_v = sin(theta_v);

		//the surface 
		theta_s = -theta_s;

		//theta_s < 0
		float cos_theta_s = cos(theta_s);
		float sin_theta_s = sin(theta_s);
		float tmp;

		for (int bounce_it = 0; bounce_it < bounce_n; bounce_it++) {

			// compute virtual surface normal
			vec3 w_face(cos_phi_sp * -sin_theta_s, sin_phi_sp * -sin_theta_s, cos_theta_s);

			// multiply the Fresnel effect at the kth intersection.
			F *= fresnelConductorExact(absDot(wi, w_face), eta, k);

			// rotate s (in counter-clockwise manner)
			tmp = sin_theta_s * cos_theta_v + cos_theta_s * sin_theta_v;
			cos_theta_s = -sin_theta_s * sin_theta_v + cos_theta_s * cos_theta_v;
			sin_theta_s = tmp;
			// - sin theta_s'   [ cos theta v -sin theta_v] [ - sin theta_s ]
			// cos theta_s' = [ sin theta v cos theta v] [ cos theta_s ]
		}

		return F;
	}
	else {

		//required variable
		//cos_phi_sp
		//sin_phi_sp
		//theta_s
		//s0
		Spectrum F(1.f);

		float theta_v = M_PI - 2 * abs(theta_s);

		float cos_theta_v = cos(theta_v);
		float sin_theta_v = sin(theta_v);

		//theta_s > 0
		float cos_theta_s = cos(theta_s);
		float sin_theta_s = sin(theta_s); // 
		float tmp;

		for (int bounce_it = 0; bounce_it < bounce_n; bounce_it++) {

			//compute virtual surface normal
			vec3 w_face(cos_phi_sp * -sin_theta_s, sin_phi_sp * -sin_theta_s, cos_theta_s);

			//multiply the Fresnel effect at the kth intersection.
			F *= fresnelConductorExact(absDot(wi, w_face), eta, k);

			// rotate s (in colockwise manner)
			tmp = sin_theta_s * cos_theta_v - cos_theta_s * sin_theta_v;
			cos_theta_s = sin_theta_s * sin_theta_v + cos_theta_s * cos_theta_v;
			sin_theta_s = tmp;
			// - sin theta_s'   [ cos -theta v -sin -theta_v] [ - sin theta_s ]
			// cos theta_s' = [ sin -theta v cos -theta v] [ cos theta_s ]
			// - sin theta_s'   [ cos theta v sin theta_v] [ - sin theta_s ]
			// cos theta_s' = [ -sin theta v cos theta v] [ cos theta_s ]
		}

		return F;
	}

}

class RoughConductorMS : public BSDF {
public:
  RoughConductorMS(const Properties &props) : BSDF(props) {
    ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

    m_specularReflectance = new ConstantSpectrumTexture(
      props.getSpectrum("specularReflectance", Spectrum(1.0f)));

    std::string materialName = props.getString("material", "Cu");

    Spectrum intEta, intK;
    if (boost::to_lower_copy(materialName) == "none") {
      intEta = Spectrum(0.0f);
      intK = Spectrum(1.0f);
    }
    else {
      intEta.fromContinuousSpectrum(InterpolatedSpectrum(
        fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
      intK.fromContinuousSpectrum(InterpolatedSpectrum(
        fResolver->resolve("data/ior/" + materialName + ".k.spd")));
    }

    Float extEta = lookupIOR(props, "extEta", "air");

    m_eta = props.getSpectrum("eta", intEta) / extEta;
    m_k = props.getSpectrum("k", intK) / extEta;

    //roughness
    MicrofacetDistribution distr(props);
    m_type = distr.getType();
    //m_sampleVisible = distr.getSampleVisible();

    m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
    if (distr.getAlphaU() == distr.getAlphaV())
      m_alphaV = m_alphaU;
    else
      m_alphaV = new ConstantFloatTexture(distr.getAlphaV());

    // scattering order
	m_scatteringOrderMax = props.getInteger("scatteringOrderMax", 10);
	m_scatteringOrderMaxForPDF = props.getInteger("scatteringOrderMaxForPDF", m_scatteringOrderMax);
    m_scatteringOrderMin = props.getInteger("scatteringOrderMin", 1);
	m_scatteringOrderMax = std::max(m_scatteringOrderMin, m_scatteringOrderMax);

  }

  RoughConductorMS(Stream *stream, InstanceManager *manager)
    : BSDF(stream, manager) {
    m_type = (MicrofacetDistribution::EType) stream->readUInt();
    m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
    m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
    m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
    m_eta = Spectrum(stream);
    m_k = Spectrum(stream);
    m_scatteringOrderMin = stream->readInt();
    m_scatteringOrderMax = stream->readInt();
    configure();
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    BSDF::serialize(stream, manager);
    stream->writeUInt((uint32_t)m_type);
    manager->serialize(stream, m_alphaU.get());
    manager->serialize(stream, m_alphaV.get());
    manager->serialize(stream, m_specularReflectance.get());
    m_eta.serialize(stream);
    m_k.serialize(stream);
    stream->writeInt(m_scatteringOrderMin);
    stream->writeInt(m_scatteringOrderMax);
  }

  void configure() {
    unsigned int extraFlags = 0;
    if (m_alphaU != m_alphaV)
      extraFlags |= EAnisotropic;

    if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
      !m_specularReflectance->isConstant())
      extraFlags |= ESpatiallyVarying;

    m_components.clear();
    m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

    /* Verify the input parameters and fix them if necessary */
    m_specularReflectance = ensureEnergyConservation(
      m_specularReflectance, "specularReflectance", 1.0f);

    m_usesRayDifferentials =
      m_alphaU->usesRayDifferentials() ||
      m_alphaV->usesRayDifferentials() ||
      m_specularReflectance->usesRayDifferentials();

    BSDF::configure();
  }

  //Compute BRDF values at given outgoing direction w_o
  Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {

    /* Stop if this component was not requested */
    if (measure != ESolidAngle ||
      Frame::cosTheta(bRec.wi) <= 0 ||
      Frame::cosTheta(bRec.wo) <= 0 ||
      ((bRec.component != -1 && bRec.component != 0) ||
        !(bRec.typeMask & EGlossyReflection)))
      return Spectrum(0.0f);

    vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
    vec3 wo(bRec.wo.x, bRec.wo.y, bRec.wo.z);

    /* Calculate the reflection half-vector */
    Vector wh = normalize(bRec.wo + bRec.wi);

    const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
    const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

    MicrofacetDistribution distr(
      m_type,
      alpha_x,
      alpha_y,
      false
    );

    float phi_sp;
    vec3 wsp, wu;
    float theta_i, theta_o;
    float abs_theta_i;
    float theta_s;
	float abs_theta_s;
    float cos_phi_sp, sin_phi_sp;

    //compute s_perpendicular (sp)
    wsp = wh;
    wsp.z = 0;
    if (wsp.length() < 0.0001f)
      wsp = vec3(1.0f, 0, 0.f);

    wsp = normalize(wsp);
    wu = cross(wsp, vec3(0, 0, 1.f));
    wu = normalize(wu);

    //compute phi of sp
    //phi_sp = atan2f(wsp.y, wsp.x);

    cos_phi_sp = wsp.x;
    sin_phi_sp = wsp.y;

    // compute thetai
    float cosbtwiwu = dot(bRec.wi, wu);
    float sinbtwiwu = sqrt(1 - cosbtwiwu*cosbtwiwu);
    theta_i = atan2f(-dot(wsp, bRec.wi), bRec.wi.z);
    theta_o = atan2f(-dot(wsp, bRec.wo), bRec.wo.z);

    float signv = 1.f;
    float pdf = 0.f;
    float D,G = 0.f;
    float theta_or, theta_ir;
    float dws_dwo;
	
    Spectrum totalReflectance(0.f), F(0.f);

//	return totalReflectance;
    theta_or = -theta_o;
    theta_ir = -theta_i; 

    //for every reflectance, we cumulate all contributions.
	for (int k = m_scatteringOrderMin; k <= m_scatteringOrderMax; k++) {

		//to compute right
		if (k & 1)
			signv = -1.f;
		else signv = 1.f;

		// case 1. a ray hits the left surface
		abs_theta_s = findSurfaceAngle(theta_o, theta_i, k, signv, true);
		float theta_v = M_PI - 2 * abs_theta_s;
		if (abs_theta_s > 0 && theta_v > 0.f && (k - 1) * theta_v < M_PI)
		{

			vec3 s0(cos_phi_sp * sin(abs_theta_s), sin_phi_sp * sin(abs_theta_s), cos(abs_theta_s));// left ( -- = + )

			//compute D
			D = distr.eval(s0);
			//compute G

			//Log(EInfo, "h: %s", wh.toString().c_str());
			computeG_k(theta_i, M_PI - abs_theta_s * 2, G, pdf , k);
			//compute F
			F = computeMultipleFresnel(cos_phi_sp, sin_phi_sp, abs_theta_s, bRec.wi, m_eta, m_k, k, true);
			//compute J
			dws_dwo = sin(abs_theta_s) / (4 * k *  absDot(bRec.wi, wh) * sqrt(1 - wh.z * wh.z));

			Assert(G >= 0.f);
			Assert(pdf >= 0.f);

			if (k == 1)
				totalReflectance += D * G * F / (4 * Frame::cosTheta(bRec.wi));
			else if (std::isfinite(dws_dwo)) totalReflectance += D * G * F * dws_dwo * absDot(s0, wi) / Frame::cosTheta(bRec.wi);

		}

		//case 2. a ray hits the right surface
		abs_theta_s = findSurfaceAngle(theta_or, theta_ir, k, signv);
		theta_v = M_PI - 2 * abs_theta_s;

		if (abs_theta_s > 0 && theta_v > 0.f && (k - 1) * theta_v < M_PI)
		{

			vec3 s0(cos_phi_sp * -sin(abs_theta_s), sin_phi_sp * -sin(abs_theta_s), cos(abs_theta_s));//right surface

			//compute D
			D = distr.eval(s0);
			//compute G
			computeG_k(theta_ir, M_PI - abs_theta_s * 2, G, pdf, k);
			//compute F
			F = computeMultipleFresnel(cos_phi_sp, sin_phi_sp, abs_theta_s, bRec.wi, m_eta, m_k, k, false);//should be changed.
			//compute J
			dws_dwo = sin(abs_theta_s) / (4 * k *  absDot(bRec.wi, wh) * sqrt(1 - wh.z * wh.z));

			Assert(G >= 0.f);
			Assert(pdf >= 0.f);

			if (k == 1)
				totalReflectance += D * G * F / (4 * Frame::cosTheta(bRec.wi));
			else if (std::isfinite(dws_dwo)) totalReflectance += D * G * F * dws_dwo * absDot(s0, wi) / Frame::cosTheta(bRec.wi);

		}
	}

	return totalReflectance;
  }

  Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {

    if (measure != ESolidAngle ||
      Frame::cosTheta(bRec.wi) < 0 ||
      Frame::cosTheta(bRec.wo) < 0 ||
      ((bRec.component != -1 && bRec.component != 0) ||
        !(bRec.typeMask & EGlossyReflection)))
      return 0.0f;

    /* Calculate the reflection half-vector */
    Vector wh = normalize(bRec.wo + bRec.wi);

    MicrofacetDistribution distr(
      m_type,
      m_alphaU->eval(bRec.its).average(),
      m_alphaV->eval(bRec.its).average(),
      false
    );
	
	float phi_sp;
	vec3 wsp, wu;
	float theta_i, theta_o;
	float abs_theta_i;
	float theta_s;
	float abs_theta_s;

	float cos_phi_sp, sin_phi_sp;

    // compute thetai
    float cosbtwiwu = dot(bRec.wi, wu);
    float sinbtwiwu = sqrt(1 - cosbtwiwu*cosbtwiwu);
    theta_i = atan2f(-dot(wsp, bRec.wi), bRec.wi.z);
    theta_o = atan2f(-dot(wsp, bRec.wo), bRec.wo.z);

    float signv = 1.f;
    float pdf = 0.f;
    float G = 0.f;
    float theta_or, theta_ir;
    float dws_dwo;
    float ratio = 0.f;
	float kpdf = 0.f;

    Spectrum totalReflectance(0.f), F(0.f);
	
    //compute s_perpendicular (sp)
    wsp = wh;
    wsp.z = 0;
    if (wsp.length() < 0.0001f)
      wsp = vec3(1.0f, 0, 0.f);

	wsp = normalize(wsp);
    wu = cross(wsp, vec3(0, 0, 1.f));
    wu = normalize(wu);

    //compute phi of sp
    //phi_sp = atan2f(wsp.y, wsp.x);

    cos_phi_sp = wsp.x;
    sin_phi_sp = wsp.y;

	theta_i = atan2f(-dot(wsp, bRec.wi), bRec.wi.z);
	theta_o = atan2f(-dot(wsp, bRec.wo), bRec.wo.z);
	
    theta_or = -theta_o;
    theta_ir = -theta_i;

    //for every reflectance, we cumulate all contributions.
    for (int k = m_scatteringOrderMin; k <= m_scatteringOrderMaxForPDF; k++) {
//	for (int k = m_scatteringOrderMin + 1; k <3; k++) {

      if (k & 1)
        signv = -1.f;
      else signv = 1.f;

      // case 1. a ray hits the left surface
	  abs_theta_s = findSurfaceAngle(theta_o, theta_i, k, signv, true);	  //find angle
	  float theta_v = M_PI - 2* abs_theta_s;

      //if (isValid(theta_o, theta_i, theta_s, k, signv)) 
	  if (abs_theta_s>0 && theta_v > 0.f && (k-1) * theta_v < M_PI)
	  {// check whether theta_s is valid or not

        vec3 s0(cos_phi_sp * sin(abs_theta_s), sin_phi_sp * sin(abs_theta_s), cos(abs_theta_s)); //TODO: check it

        //compute ratio
		computeG_k(theta_i, M_PI - abs_theta_s * 2, G, ratio, k);

		//ratio = abs(ratio);

        //compute J
        if ( k== 1 )
          dws_dwo = 1 / ( 4.f *  absDot( bRec.wi, wh ) );
        else dws_dwo = sin(abs_theta_s) / (4.f * k *  absDot( bRec.wi, wh ) * sqrt( 1 - wh.z * wh.z ));

		kpdf = distr.pdf(bRec.wi, s0) * dws_dwo * ratio;

		Assert(G >= 0.f);
		Assert(ratio >= 0.f);

		if (std::isfinite(kpdf)) {
			pdf += kpdf;
		}

      }

      //case 2. a ray hits the right surface
	  abs_theta_s = findSurfaceAngle(theta_or, theta_ir, k, signv);
	  theta_v = M_PI - 2 * abs_theta_s;

     // if (isValid(theta_or, theta_ir, theta_s, k, signv)) 
	  if(abs_theta_s>0 && theta_v > 0.f && (k - 1) * theta_v < M_PI)
	  {

        vec3 s0(cos_phi_sp * sin(-abs_theta_s), sin_phi_sp * sin(-abs_theta_s), cos(abs_theta_s)); //the direction should be filpped but doesnt matter to compute the ratio.

        //compute ratio
		computeG_k(theta_ir, M_PI - abs_theta_s * 2, G, ratio, k);

		//ratio = abs(ratio);
		//compute J
        if (k == 1)
          dws_dwo = 1 / (4.f *  absDot(bRec.wi, wh));
        else dws_dwo = sin(abs_theta_s) / (4.f * k *  absDot(bRec.wi, wh) * sqrt(1 - wh.z * wh.z));

		kpdf = distr.pdf(bRec.wi, s0) * dws_dwo * ratio;

		Assert(G >= 0.f);
		Assert(ratio >= 0.f);

		if (std::isfinite(kpdf)) {
			pdf += kpdf;
		}

      }
    }

    return pdf;
  }
  
  Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {

	  if (Frame::cosTheta(bRec.wi) < 0 ||
		  ((bRec.component != -1 && bRec.component != 0) ||
			  !(bRec.typeMask & EGlossyReflection)))
		  return Spectrum(0.0f);

	  Float pdf;
	  return this->sample(bRec, pdf, sample);

  }

  Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
	  if (Frame::cosTheta(bRec.wi) < 0 ||
		  ((bRec.component != -1 && bRec.component != 0) ||
			  !(bRec.typeMask & EGlossyReflection)))
		  return Spectrum(0.0f);

	  /* Construct the microfacet distribution matching the
	  roughness values at the current surface position. */
	  MicrofacetDistribution distr(
		  m_type,
		  m_alphaU->eval(bRec.its).average(),
		  m_alphaV->eval(bRec.its).average(),
		  false
	  );

	  /* Sample M, the microfacet normal */
	  float microfacetPDF;
	  Normal s0 = distr.sample(bRec.wi, sample, microfacetPDF);

	  if (microfacetPDF == 0.f)
		  return Spectrum(0.0f);

	  pdf = microfacetPDF;

	  float phi_sp;
	  vec3 wsp, wu;
	  vec3 wi_projected, wo_projected;
	  float theta_i, theta_o;
	  float theta_v, half_theta_v;
	  float abs_theta_i;
	  float theta_s;

	  float cos_phi_sp, sin_phi_sp;
	  float cos_theta_v, sin_theta_v;
	  float cos_theta_s, sin_theta_s;
	  float cos_theta_s0, sin_theta_s0;

	  //compute s_perpendicular (sp)
	  wsp = s0;
	  wsp.z = 0;
	  if (wsp.length() < 0.0001f)
		  wsp = vec3(1.0f, 0, 0.f);
	  wsp = normalize(wsp);
	  wu = cross(wsp, vec3(0, 0, 1.f));
	  wu = normalize(wu);

	  cos_phi_sp = wsp.x;
	  sin_phi_sp = wsp.y;

	  // compute thetai
	  float cosbtwiwu = dot(bRec.wi, wu);
	  float sinbtwiwu = sqrt(1 - cosbtwiwu*cosbtwiwu);
	  theta_i = atan2f(-dot(wsp, bRec.wi), bRec.wi.z);


	  theta_s = -acosf(s0.z); //theta_s < 0
	  theta_v = M_PI + 2 * theta_s; // 2 * ( pi/2 - | theta_s | )
	  half_theta_v = 0.5f * theta_v;

	  cos_theta_v = cos(theta_v);
	  sin_theta_v = sin(theta_v);

	  int bounce_cnt;
	  float G;

	  //sample direction and compute G 
	  pdf = pdf * sampleDirectionAndComputeG(theta_i, theta_v, bounce_cnt, G, m_sampleVisible);
	  
	  Assert(G >= 0.f);
	  Assert(pdf >= 0.f);

	  if (m_scatteringOrderMin > bounce_cnt || bounce_cnt > m_scatteringOrderMax)
		  return Spectrum(0.f);

	  //compute F
	  Spectrum F(1.f);

	  cos_theta_s0 = s0.z;
	  sin_theta_s0 = -sqrt(1 - s0.z*s0.z); //theta_s < 0 , so sin < 0
	  vec3 ws0(cos_phi_sp * -sin_theta_s0, sin_phi_sp * -sin_theta_s0, cos_theta_s0);
	  cos_theta_s = cos_theta_s0;
	  sin_theta_s = sin_theta_s0;

	  theta_o = theta_i + M_PI - bounce_cnt * theta_v;

	  if (bounce_cnt & 1)
		  theta_o = -1 * theta_o;

	  vec3 w_op(cos_phi_sp * -sin(theta_o), sin_phi_sp * -sin(theta_o), cos(theta_o));

	  if (w_op.z <= 0.f)
		  return Spectrum(0.0f);

	  bRec.wo = normalize(sinbtwiwu * w_op - cosbtwiwu * wu);
	  bRec.eta = 1.0f;
	  bRec.sampledComponent = 0;
	 // bRec.sampledType = EDeltaReflection; //TURN OFF MIS during BSDF sampling.
	   bRec.sampledType = EGlossyReflection;
	  vec3 wh = normalize(bRec.wi + bRec.wo);

	  /* Side check */
	  if (Frame::cosTheta(bRec.wo) <= 0)
		  return Spectrum(0.0f);

	  F = computeMultipleFresnel(cos_phi_sp, sin_phi_sp, theta_s, bRec.wi, m_eta, m_k, bounce_cnt, true);// hit the left surface

	  float model;

	  float dws_dwo = abs((sin_theta_s0) / (4 * bounce_cnt *  absDot(bRec.wo, wh) * sqrt(1 - wh.z * wh.z)));

	 model = distr.eval(s0) *  G  * absDot(bRec.wi, s0) / (pdf * Frame::cosTheta(bRec.wi));
	 pdf = this->pdf(bRec, ESolidAngle);//To compute MIS weight
	 
	 if (pdf <= 0.f) {
		 pdf = 1.f;
		 return Spectrum(0.f);
	 }
	  return F * model;
  }

  Float getRoughness(const Intersection &its, int component) const {
      return 0.5f * (m_alphaU->eval(its).average()
        + m_alphaV->eval(its).average());
    }

  void addChild(const std::string &name, ConfigurableObject *child) {
      if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
        if (name == "alpha")
          m_alphaU = m_alphaV = static_cast<Texture *>(child);
        else if (name == "alphaU")
          m_alphaU = static_cast<Texture *>(child);
        else if (name == "alphaV")
          m_alphaV = static_cast<Texture *>(child);
        else if (name == "specularReflectance")
          m_specularReflectance = static_cast<Texture *>(child);
        else
          BSDF::addChild(name, child);
      }
      else {
        BSDF::addChild(name, child);
      }
    }

  std::string toString() const {
    std::ostringstream oss;
    oss << "RoughConductorMS[" << endl
      << "  id = \"" << getID() << "\"," << endl
      << "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
		<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
		<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
		<< "  scateringOrderMin = " << m_scatteringOrderMin << "," << endl
		<< "  scatteringOrderMax = " << m_scatteringOrderMax << "," << endl
		<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
      << "  eta = " << m_eta.toString() << "," << endl
      << "  k = " << m_k.toString() << endl
      << "]";
    return oss.str();
  }

  Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()

private:
	MicrofacetDistribution::EType m_type;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	Spectrum m_eta, m_k;
  int m_scatteringOrderMax;
  int m_scatteringOrderMaxForPDF;
  int m_scatteringOrderMin;
  bool m_sampleVisible;
};

/**
* GLSL port of the rough conductor shader. This version is much more
* approximate -- it only supports the Ashikhmin-Shirley distribution,
* does everything in RGB, and it uses the Schlick approximation to the
* Fresnel reflectance of conductors. When the roughness is lower than
* \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
* reasonably well in a VPL-based preview.
*/

class RoughConductorMSShader : public Shader {
public:
  RoughConductorMSShader(Renderer *renderer, const Texture *specularReflectance,
    const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
    const Spectrum &k) : Shader(renderer, EBSDFShader),
    m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV) {
    m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
    m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
    m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

    /* Compute the reflectance at perpendicular incidence */
    m_R0 = fresnelConductorExact(1.0f, eta, k);
  }

  bool isComplete() const {
    return m_specularReflectanceShader.get() != NULL &&
      m_alphaUShader.get() != NULL &&
      m_alphaVShader.get() != NULL;
  }

  void putDependencies(std::vector<Shader *> &deps) {
    deps.push_back(m_specularReflectanceShader.get());
    deps.push_back(m_alphaUShader.get());
    deps.push_back(m_alphaVShader.get());
  }

  void cleanup(Renderer *renderer) {
    renderer->unregisterShaderForResource(m_specularReflectance.get());
    renderer->unregisterShaderForResource(m_alphaU.get());
    renderer->unregisterShaderForResource(m_alphaV.get());
  }

  void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
    parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
  }

  void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
    program->setParameter(parameterIDs[0], m_R0);
  }

  void generateCode(std::ostringstream &oss,
    const std::string &evalName,
    const std::vector<std::string> &depNames) const {
    oss << "uniform vec3 " << evalName << "_R0;" << endl
      << endl
      << "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
      << "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
      << "    if (ds <= 0.0)" << endl
      << "        return 0.0f;" << endl
      << "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
      << "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
      << "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
      << "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
      << "}" << endl
      << endl
      << "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
      << "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
      << "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
      << "        return 0.0;" << endl
      << "    float nDotM = cosTheta(m);" << endl
      << "    return min(1.0, min(" << endl
      << "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
      << "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
      << "}" << endl
      << endl
      << "vec3 " << evalName << "_schlick(float ct) {" << endl
      << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
      << "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
      << "}" << endl
      << endl
      << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
      << "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
      << "    	return vec3(0.0);" << endl
      << "   vec3 H = normalize(wi + wo);" << endl
      << "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
      << "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
      << "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
      << "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
      << "   float G = " << evalName << "_G(H, wi, wo);" << endl
      << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
      << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
      << "}" << endl
      << endl
      << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
      << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
      << "    	return vec3(0.0);" << endl
      << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);" << endl
      << "}" << endl;
  }
  MTS_DECLARE_CLASS()
private:
  ref<const Texture> m_specularReflectance;
  ref<const Texture> m_alphaU;
  ref<const Texture> m_alphaV;
  ref<Shader> m_specularReflectanceShader;
  ref<Shader> m_alphaUShader;
  ref<Shader> m_alphaVShader;
  Spectrum m_R0;
};

Shader *RoughConductorMS::createShader(Renderer *renderer) const {
  return new RoughConductorMSShader(renderer,
    m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(RoughConductorMSShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughConductorMS, false, BSDF)
MTS_EXPORT_PLUGIN(RoughConductorMS, "Rough conductor MS BRDF (new)");
MTS_NAMESPACE_END
