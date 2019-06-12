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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/chisquare.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/testcase.h>
#include <boost/bind.hpp>

#define SIGNIFICANCE_LEVEL 0.0025f

#if defined(SINGLE_PRECISION)
#define ERROR_REQ 1e-2f
#else
#define ERROR_REQ 1e-5
#endif

MTS_NAMESPACE_BEGIN


class TestOurTest : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
		MTS_DECLARE_TEST(test01_EnergyConservation)
		MTS_DECLARE_TEST(test02_Reciprocity)
		MTS_DECLARE_TEST(test03_SingleScattering)
		MTS_END_TESTCASE()

		/**
		* Replayable fake sampler
		*/
		class FakeSampler : public Sampler {
		public:
			FakeSampler(Sampler *sampler)
				: Sampler(Properties()), m_sampler(sampler) { }

			Float next1D() {
				while (m_sampleIndex >= m_values.size())
					m_values.push_back(m_sampler->next1D());
				return m_values[m_sampleIndex++];
			}

			Point2 next2D() {
				return Point2(next1D(), next1D());
			}

			void clear() {
				m_values.clear();
				m_sampleIndex = 0;
			}

			void rewind() {
				m_sampleIndex = 0;
			}

			ref<Sampler> clone() {
				SLog(EError, "Not supported!");
				return NULL;
			}

			std::string toString() const { return "FakeSampler[]"; }
		private:
			ref<Sampler> m_sampler;
			std::vector<Float> m_values;
	};


	class BSDFAdapter {
	public:
		BSDFAdapter(const BSDF *bsdf, Sampler *sampler, const Vector &wi, int component)
			: m_bsdf(bsdf), m_sampler(sampler), m_wi(wi), m_component(component),
			m_largestWeight(0) {
			m_fakeSampler = new FakeSampler(m_sampler);
			m_its.uv = Point2(0.0f);
			m_its.dpdu = Vector(1, 0, 0);
			m_its.dpdv = Vector(0, 1, 0);
			m_its.dudx = m_its.dvdy = 0.01f;
			m_its.dudy = m_its.dvdx = 0.00f;
			m_its.shFrame = Frame(Normal(0, 0, 1));
			//			m_isSymmetric = true;
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			Point2 sample(m_sampler->next2D());
			BSDFSamplingRecord bRec(m_its, m_fakeSampler);
			bRec.mode = EImportance;
			bRec.component = m_component;
			bRec.wi = m_wi;

#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
#endif

			Float pdfVal, sampledPDF;

			/* Check the various sampling routines for agreement
			amongst each other */
			m_fakeSampler->clear();
			Spectrum sampled = m_bsdf->sample(bRec, sampledPDF, sample);
			m_fakeSampler->rewind();
			Spectrum sampled2 = m_bsdf->sample(bRec, sample);
			EMeasure measure = ESolidAngle;
			if (!sampled.isZero())
				measure = BSDF::getMeasure(bRec.sampledType);

			if (sampled.isZero() && sampled2.isZero())
				return boost::make_tuple(Vector(0.0f), 0.0f, measure);

			Spectrum f = m_bsdf->eval(bRec, measure);
			pdfVal = m_bsdf->pdf(bRec, measure);
			Spectrum manual = f / pdfVal;

#if 0
			if (m_isSymmetric) {
				/* Check for non-symmetry */
				BSDFSamplingRecord bRecRev(bRec);
				bRecRev.reverse();
				bRec.mode = EImportance;
				Spectrum fFwd = f;
				Spectrum fRev = m_bsdf->eval(bRecRev, measure);
				if (measure == ESolidAngle) {
					fFwd /= std::abs(Frame::cosTheta(bRec.wo));
					fRev /= std::abs(Frame::cosTheta(bRecRev.wo));
				}
				Float max = std::max(fFwd.max(), fRev.max());
				if (max > 0) {
					Float err = (fFwd - fRev).max() / max;
					if (err > Epsilon) {
						Log(EWarn, "Non-symmetry in %s: %s vs %s, %s", m_bsdf->toString().c_str(),
							fFwd.toString().c_str(), fRev.toString().c_str(), bRec.toString().c_str());
						m_isSymmetric = false;
					}
				}
			}
#endif

			if (!sampled.isValid() || !sampled2.isValid() || !manual.isValid()) {
				Log(EWarn, "Oops: sampled=%s, sampled2=%s, manual=%s, sampledPDF=%f, "
					"pdf=%f, f=%s, bRec=%s, measure=%i", sampled.toString().c_str(),
					sampled2.toString().c_str(), manual.toString().c_str(),
					sampledPDF, pdfVal, f.toString().c_str(), bRec.toString().c_str(),
					measure);
				return boost::make_tuple(bRec.wo, 0.0f, ESolidAngle);
			}

			bool mismatch = false;
			for (int i = 0; i < SPECTRUM_SAMPLES; ++i) {
				Float a = sampled[i], b = sampled2[i], c = manual[i];
				Float min = std::min(std::min(a, b), c);
				Float err = std::max(std::max(std::abs(a - b), std::abs(a - c)), std::abs(b - c));
				m_largestWeight = std::max(m_largestWeight, a);

				if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
					mismatch = true;
				else if (min > ERROR_REQ && err / min > ERROR_REQ) // relative error threshold
					mismatch = true;
			}

			if (mismatch)
				Log(EWarn, "Potential inconsistency: sampled=%s, sampled2=%s, manual=%s, sampledPDF=%f, "
					"pdf=%f, f=%s, bRec=%s, measure=%i", sampled.toString().c_str(),
					sampled2.toString().c_str(), manual.toString().c_str(),
					sampledPDF, pdfVal, f.toString().c_str(), bRec.toString().c_str(),
					measure);

			mismatch = false;
			Float min = std::min(pdfVal, sampledPDF);
			Float err = std::abs(pdfVal - sampledPDF);

			if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
				mismatch = true;
			else if (min > ERROR_REQ && err / min > ERROR_REQ) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Potential inconsistency: pdfVal=%f, sampledPDF=%f",
					pdfVal, sampledPDF);

#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
#endif

			return boost::make_tuple(bRec.wo, 1.0f, measure);
		}

		Spectrum eval(const Vector &wo, EMeasure measure) {
			BSDFSamplingRecord bRec(m_its, m_wi, wo);
			bRec.mode = ERadiance;
			bRec.component = m_component;

#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
#endif

			Spectrum f = m_bsdf->eval(bRec, measure);

			if (f.isZero())
				return Spectrum(0.0f);

#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
#endif

			return f;
		}

		Float pdf(const Vector &wo, EMeasure measure) {
			BSDFSamplingRecord bRec(m_its, m_wi, wo);
			bRec.mode = EImportance;
			bRec.component = m_component;

#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
#endif

			if (m_bsdf->eval(bRec, measure).isZero())
				return 0.0f;

			Float result = m_bsdf->pdf(bRec, measure);

#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
#endif
			return result;
		}

		inline Float getLargestWeight() const { return m_largestWeight; }
	private:
		Intersection m_its;
		ref<const BSDF> m_bsdf;
		ref<Sampler> m_sampler;
		ref<FakeSampler> m_fakeSampler;
		Vector m_wi;
		int m_component;
		Float m_largestWeight;
	};

	class EmitterAdapter {
	public:
		EmitterAdapter(const Emitter *emitter, Sampler *sampler)
			: m_emitter(emitter), m_sampler(sampler), m_pRec(0.0f) {
			emitter->samplePosition(m_pRec, m_sampler->next2D());
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
#endif

			DirectSamplingRecord dRec(Point(0.0f), 0);
			m_emitter->sampleDirect(dRec, m_sampler->next2D());

#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
#endif

			return boost::make_tuple(dRec.d, 1.0f, dRec.measure);
		}

		Float pdf(const Vector &d, EMeasure measure) const {
			if (measure != ESolidAngle)
				return 0.0f;

			DirectSamplingRecord dRec(Point(0.0f), 0);
			dRec.d = d;
			dRec.measure = ESolidAngle;

#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
#endif

			Float result = m_emitter->pdfDirect(dRec);

#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
#endif

			return result;
		}

	private:
		ref<const Emitter> m_emitter;
		ref<Sampler> m_sampler;
		PositionSamplingRecord m_pRec;
	};

	void test01_EnergyConservation() {

		/* Load a set of BSDF instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_bsdfours.xml");
		ref<Scene> scene = loadScene(scenePath);

		const ref_vector<ConfigurableObject> &objects = scene->getReferencedObjects();
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		Log(EInfo, "Verifying BSDF sampling routines ..");

		int wiSamples = 20;
		int i;
		i = 0;
		if (objects[i]->getClass()->derivesFrom(MTS_CLASS(BSDF))) {

			const BSDF *bsdf = static_cast<const BSDF *>(objects[i].get());

			for (int j = 0; j < wiSamples; ++j) {
				Vector wi;

				if (bsdf->getType() & BSDF::EBackSide)
					wi = warp::squareToUniformSphere(sampler->next2D());
				else
					wi = warp::squareToCosineHemisphere(sampler->next2D());

				BSDFAdapter adapter(bsdf, sampler, wi, -1);

				Spectrum value_quadrature(0.f);
				for (double theta_o = 0; theta_o < M_PI; theta_o += 0.005)
					for (double phi_o = 0; phi_o < 2.0*M_PI; phi_o += 0.005)
					{
						const Vector wo(cos(phi_o)*sin(theta_o), sin(phi_o)*sin(theta_o), cos(theta_o));
						// stochastic evaluation
						const int N = 1;
						Spectrum value_current(0.f);
						for (int n = 0; n < N; ++n)
						{
							value_current += adapter.eval(wo, ESolidAngle) / (double)N;
						}
						value_quadrature += 0.005*0.005*abs(sin(theta_o)) * value_current;
					}

				// display
				Log(EInfo, "int f_p(wi = %s, wo ) dwo = \t\t %f", wi.toString().c_str(), value_quadrature[0]);
			}
		}
	}

	void test02_Reciprocity() {


		/* Load a set of BSDF instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_bsdfours.xml");
		ref<Scene> scene = loadScene(scenePath);

		const ref_vector<ConfigurableObject> &objects = scene->getReferencedObjects();
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		Log(EInfo, "Verifying BSDF sampling routines ..");

		int wiSamples = 100;
		int i;
		i = 0;
		if (objects[i]->getClass()->derivesFrom(MTS_CLASS(BSDF))) {

			const BSDF *bsdf = static_cast<const BSDF *>(objects[i].get());

			float totalerr = 0.f;
			float err = 0.f;

			for (int j = 0; j < wiSamples; ++j) {
				Vector wi;
				Vector wo;

				if (bsdf->getType() & BSDF::EBackSide) {
					wi = warp::squareToUniformSphere(sampler->next2D());
					wo = warp::squareToUniformSphere(sampler->next2D());
				}
				else {
					wi = warp::squareToCosineHemisphere(sampler->next2D());
					wo = warp::squareToCosineHemisphere(sampler->next2D());
				}

				BSDFAdapter adapteri(bsdf, sampler, wi, -1);
				BSDFAdapter adaptero(bsdf, sampler, wo, -1);

				Spectrum valueio = adapteri.eval(wo, ESolidAngle)* wi.z;
				Spectrum valueoi = adaptero.eval(wi, ESolidAngle)* wo.z;

				err = 0.f;
				for (int si = 0; si < SPECTRUM_SAMPLES; ++si) {
					err += (valueio[si] - valueoi[si])*(valueio[si] - valueoi[si]);
				}
				Log(EInfo, "%s %s Value_io, Value_oi:\t\t%s, %s\terror:\t\t%f", wi.toString().c_str(), wo.toString().c_str(), valueio.toString().c_str(), valueoi.toString().c_str(), err);
				totalerr += err;
			}
			Log(EInfo, "Average error:\t\t%f", totalerr / wiSamples);
		}
	}

	void test03_SingleScattering() {

		/* Load a set of BSDF instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_bsdfours.xml");
		ref<Scene> scene = loadScene(scenePath);

		const ref_vector<ConfigurableObject> &objects = scene->getReferencedObjects();
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		Log(EInfo, "Verifying BSDF sampling routines ..");

		int wiSamples = 100;

		if (objects[1]->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (objects[2]->getClass()->derivesFrom(MTS_CLASS(BSDF))) {



				const BSDF *bsdf0 = static_cast<const BSDF *>(objects[1].get());// our model in single reflection
				const BSDF *bsdf1 = static_cast<const BSDF *>(objects[2].get());// V-groove model

				Log(EInfo, "Two BSDF name:\t\t%s, %s", bsdf0->toString().c_str(), bsdf1->toString().c_str());
				
				float totalerr = 0.f;
				float err = 0.f;

				for (int j = 0; j < wiSamples; ++j) {
					Vector wi;
					Vector wo;

					wi = warp::squareToCosineHemisphere(sampler->next2D());
					wo = warp::squareToCosineHemisphere(sampler->next2D());

					BSDFAdapter adapter0(bsdf0, sampler, wi, -1);
					BSDFAdapter adapter1(bsdf1, sampler, wi, -1);

					Spectrum value0 = adapter0.eval(wo, ESolidAngle);
					Spectrum value1 = adapter1.eval(wo, ESolidAngle);

					Spectrum diff = value0 - value1;
					
					err = 0.f;
					for (int si = 0; si < SPECTRUM_SAMPLES; ++si) {
						err += (value0[si] - value1[si])*(value0[si] - value1[si]);
					}

					Log(EInfo, "%s %s Value_ourSingle, Value_Single:\t\t%s, %s\terror:\t\t%f", wi.toString().c_str(), wo.toString().c_str(), value0.toString().c_str(), value1.toString().c_str(),err);
					totalerr += err;
				}
				Log(EInfo, "Average error:\t\t%f", totalerr / wiSamples);

			}
		}
	}
};

MTS_EXPORT_TESTCASE(TestOurTest, "Our BSDF test (Energy conservation, reciprocity and single scattering)")
MTS_NAMESPACE_END
