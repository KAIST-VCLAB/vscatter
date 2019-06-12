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

#if !defined(__VSCATTER_H__)
#define __VSCATTER_H__

#define vec3 Vector
#define vec2 Vector

///////////////////////////////////////////////////////////////////////////////

static inline float generateRandomNumber(){

  const float U = ((float)rand()) / (float)RAND_MAX;
	return U;

}

static inline float sampleDirectionAndComputeG(float &theta_i, float &theta_v, int &bounce_cnt, float &G, const bool &Visible) {

  int valid_bounce_cnt;
  float li_x, li_y;
  float fa, fb, fc;
  float half_theta_v;
  float theta_tmp;
  float pdf;

  half_theta_v = 0.5f * theta_v;

  li_x = cos( theta_i );
  li_y = sin( theta_i );

  if (theta_i <= -half_theta_v) { //hit only left side

    valid_bounce_cnt = floorf((M_PI + 2.f * theta_i) / theta_v) + 1;

    fa = li_x * sin(-1 * half_theta_v) + li_y * cos(half_theta_v);
    fb = li_x * sin(half_theta_v) + li_y * cos(half_theta_v);

    //k th edge point
    theta_tmp = ((valid_bounce_cnt - 1.f)  * theta_v + half_theta_v);
    fc = li_x * -sin(theta_tmp) + li_y * cos(theta_tmp);


    const float U1 = (fb - fa) * generateRandomNumber();

    if (U1 < fc - fa) {
      bounce_cnt = valid_bounce_cnt-1;
      G = (fa - fc) / fa;
      pdf = (fc - fa) / (fb - fa);
    }

    else {
      bounce_cnt = valid_bounce_cnt;
      G = (fc - fb) / fa;
      pdf = (fb - fc) / (fb - fa);
    }
  }
  
  else if (theta_i < half_theta_v) {

    valid_bounce_cnt = floorf((M_PI + theta_i) / theta_v + 0.5f);
    //fa<0
    fa = li_x * -sin(half_theta_v) + li_y * cos(half_theta_v);
    theta_tmp = ((valid_bounce_cnt - 1) * theta_v + half_theta_v);
    fc = li_x * -sin(theta_tmp) + li_y * cos(theta_tmp);

    const float U1 = fa * generateRandomNumber();

     if (0 > fc && fc > fa) {
       if (U1 > fc) {
         bounce_cnt = valid_bounce_cnt;
         G = fc / fa;
         pdf = G;

        }
       else{
         bounce_cnt = valid_bounce_cnt - 1;
         G = (fa - fc) / fa;
         pdf = G;
       }
    }
    else if(fc >= 0){
      bounce_cnt = valid_bounce_cnt - 1;
      G = 1;
      pdf = 1;
    }
    else {
      bounce_cnt = valid_bounce_cnt;
      G = 1;
      pdf = 1;
    }
  }

  else {
    G = 0.;
    pdf = 1.;
  }

  return pdf;

}

///////////////////////////////////////////////////////////////////////////////

static inline float signf(float x) {

	return (float(0) < x) - (x < float(0));

}

static inline void computeG_k(float theta_i, float theta_v, float &G, float &pdf,  int bounce_cnt) {

  int valid_bounce_cnt;
  float li_x, li_y;
  float fa, fb, fc;
  float half_theta_v;
  float theta_tmp;
  float coverage;

  half_theta_v = 0.5f * theta_v;

  G = 0.f;
  pdf = 1.f;

  li_x = cos(theta_i);
  li_y = sin(theta_i);

  if (theta_i <= -half_theta_v) { //hit only left side

    valid_bounce_cnt = floorf((M_PI + 2.f * theta_i) / theta_v) + 1;

    fa = li_x * sin(-1 * half_theta_v) + li_y * cos(half_theta_v);
    fb = li_x * sin(half_theta_v) + li_y * cos(half_theta_v);
	coverage = fa - fb;

    //k th edge point
    theta_tmp = ((valid_bounce_cnt - 1.f)  * theta_v + half_theta_v);
    fc = li_x * -sin(theta_tmp) + li_y * cos(theta_tmp);
	
	if ( fa<=fc && fc<=fb) {
		if (bounce_cnt == valid_bounce_cnt - 1) {
			G = (fa - fc) / fa;
			pdf = (fa - fc) / coverage;
		}
		else if (bounce_cnt == valid_bounce_cnt) {
			G = (fc - fb) / fa;
			pdf = (fc - fb) / coverage;
		}
	}
	else if (fc > fb) {
		if (bounce_cnt == valid_bounce_cnt - 1) {
			pdf = 1.f;
			G = (fa - fb) / fa;

		}
	}
	else {
		if (bounce_cnt == valid_bounce_cnt) {
			pdf = 1.f;
			G = (fa - fb) / fa;
		}
	}
  }

  else if (theta_i < half_theta_v) {// hit both surfaces 

    valid_bounce_cnt = floorf((M_PI + theta_i) / theta_v + 0.5f);
    //fa<0
    fa = li_x * -sin(half_theta_v) + li_y * cos(half_theta_v);
    theta_tmp = ((valid_bounce_cnt - 1) * theta_v + half_theta_v);
    fc = li_x * -sin(theta_tmp) + li_y * cos(theta_tmp);
	coverage = fa;

    if (0 > fc && fc > fa) {
		if (bounce_cnt == valid_bounce_cnt) {
			pdf = G = fc / fa;
		}
		else if (bounce_cnt == valid_bounce_cnt - 1) {
			pdf = G =  (fa - fc) / fa;
		}
    }
    else if (fc >= 0) {
		if (bounce_cnt == valid_bounce_cnt - 1) {
			pdf = G =  1.f;
		}
    }
    else {
      if (bounce_cnt == valid_bounce_cnt)
		  pdf = G = 1.f;
	}
  }
}

static inline float findSurfaceAngle(float &theta_o, float &theta_i, int bounce_cnt, float signv, bool isleft = true) {

  if ( isleft )
	  //return 0.5*((theta_i + M_PI - signv * theta_o) / bounce_cnt - M_PI);
	      return 0.5*((signv * theta_o - theta_i - M_PI) / bounce_cnt + M_PI);
  //else return 0.5*((-signv * theta_o + theta_i - M_PI) / bounce_cnt + M_PI);

}

#endif
