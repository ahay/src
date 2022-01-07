/* 3D Hybrid Radon transform - diffractions + reflections in the time dip-angle domain */
/*
  Copyright (C) 2013 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#include "stretch.h"

// data parameters
static int tn_;    static float to_, td_;
static int dipn_;  static float dipo_, dipd_;
static int sdipn_; static float sdipo_, sdipd_;
// model parameters
static int xin_;    static float xio_, xid_;
static int sxin_;   static float sxio_, sxid_;
static int dip0n_;  static float dip0o_, dip0d_;
static int sdip0n_; static float sdip0o_, sdip0d_;

static float *tmp, *amp, *str, *tx;

static float *tableRefl_, *tableDiff_;

static bool isAA_;
static int invMod_;
static int tLim_; 

void ditime3d_init (float dipo,   float dipd,   int dipn,   // x-dip angle axis 
					float sdipo,  float sdipd,  int sdipn,  // y-dip angle axis 
				    float xio,    float xid,    int xin,    // xi axis in x-direction
				    float sxio,   float sxid,   int sxin,   // xi axis in y-direction
				    float dip0o,  float dip0d,  int dip0n,  // x-refl dip axis
				    float sdip0o, float sdip0d, int sdip0n, // y-refl dip axis
				    float to,     float td,     int tn,     // time axis 
				    bool isAA,                              // antialiasing
				    int invMod) 						 
/*< initialize >*/
{
    int reflSize, diffSize;
    float CONVPARAM;
    float *pTableR, *pTableD;
    int id, ixi, id0, it, sid, sixi, sid0;
    float curDip, a, tan_a, tan_sa, curXi, scurXi, aux_diff, curTime, a0, sa0, tan_a0=0, tan_sa0, aux_refl, sa;
	float curDip0, scurDip, scurDip0;

	float k2, ke, k0, k1, kM;

    to_   = to;   td_   = td;   tn_   = tn;  

    dipo_ = dipo; dipd_ = dipd; dipn_ = dipn;  
    sdipo_ = sdipo; sdipd_ = sdipd; sdipn_ = sdipn;  

    dip0o_ = dip0o; dip0d_ = dip0d; dip0n_ = dip0n;  
    sdip0o_ = sdip0o; sdip0d_ = sdip0d; sdip0n_ = sdip0n;  
    xio_  = xio; xid_  = xid;  xin_  = xin;  
    sxio_ = sxio; sxid_ = sxid; sxin_ = sxin;  

    tLim_ = to_ + td_ * (tn_ - 1);

    invMod_ = invMod;
    isAA_   = isAA;

    if (isAA_) {
		sf_aastretch_init  (false, tn_, to_, td_, tn_);
		sf_halfint_init (true, 2 * tn_, 1.f - 1.f / tn_);
    } else {
		stretch_init  (tn_, to_, td_, tn_);
	}

	amp = sf_floatalloc (tn_);
    str = sf_floatalloc (tn_);
    tx  = sf_floatalloc (tn_);
    tmp = sf_floatalloc (tn_);

    // shift tables
    reflSize = tn_ * dip0n_ * sdip0n_ * dipn_ * sdipn_;
    diffSize = tn_ * xin_ * dipn_ * sxin_ * sdipn_;

    tableRefl_ = sf_floatalloc (reflSize);
    tableDiff_ = sf_floatalloc (diffSize);
	
    CONVPARAM = SF_PI / 180.f;

    pTableR = tableRefl_;
    pTableD = tableDiff_;

    for (sid = 0; sid < sdipn_; ++sid) { 
		scurDip = sdipo_ + sid * sdipd_;
		sa      = scurDip * CONVPARAM;
		tan_sa  = tan (sa);
	    for (id = 0; id < dipn_; ++id) { 
			curDip = dipo_ + id * dipd_;
			a      = curDip * CONVPARAM;
			tan_a  = tan (a);
	
			// diffraction part
			for (sixi = 0; sixi < sxin_; ++sixi) { 
			    scurXi = sxio_ + sixi * sxid_;
				for (ixi = 0; ixi < xin_; ++ixi) { 
				    curXi = xio_ + ixi * xid_;

					k2 = scurXi * tan_sa + curXi * tan_a;
					ke = scurXi*scurXi + curXi*curXi + 1.f;
	
				    aux_diff = k2 + sqrt (k2 * k2 + ke);
	
				    for (it = 0; it < tn_; ++it, ++pTableD) { 
						curTime = to_ + it * td_;
						*pTableD = curTime * aux_diff;			    
				    }
				}
			}		
				
			// reflection part
			for (sid0 = 0; sid0 < sdip0n_; ++sid0) { 
			    scurDip0 = sdip0o_ + sid0 * sdip0d_;
			    sa0 = scurDip0 * CONVPARAM;
			    tan_sa0 = tan (sa0);	

				k0 = sqrt (1.f + tan_a0*tan_a0 + tan_sa0*tan_sa0);

				for (id0 = 0; id0 < dip0n_; ++id0) { 
				    curDip0 = dip0o_ + id0 * dip0d_;
				    a0 = curDip0 * CONVPARAM;
				    tan_a0 = tan (a0);	
	
					kM = sqrt (1.f + tan_a*tan_a + tan_sa*tan_sa);
					k1 = tan_sa0*tan_sa + tan_a0*tan_a;
			
				    aux_refl = 1.f / (k0 * kM - k1);
		
				    for (it = 0; it < tn_; ++it, ++pTableR) {		
						curTime = to_ + it * td_;
						*pTableR = curTime * aux_refl;		
				    }
				}	
	    	}	

		}
	}

    return;
}

void ditime3d_close (void)
/*< free allocated storage >*/
{
    free(amp);
    free(str);
    free(tx);
    free(tmp);

    if (isAA_) {
		sf_aastretch_close();
		sf_halfint_close();
    } else {
		stretch_close();
	}
	
    free (tableRefl_);
    free (tableDiff_);

    return;
}

void ditime3d_lop (bool adj, bool add, int modelSize, int dataSize, 
				   float* modl, float* data) 
/*< operator >*/
{
    float *pTableR, *pTableD;
    int id, ixi, it, sid, sixi, offset, id0, sid0;
    float curTime, t;

    switch (invMod_) {
	case 0:
	    if (modelSize != tn_ * xin_ * sxin_ || dataSize != tn_ * dipn_ * sdipn_) 
		sf_error ("%s: wrong dimensions modelSize=%d dataSize=%d",__FILE__, modelSize, dataSize);
	    break;
	case 1:
	    if (modelSize != tn_ * dip0n_ * sdip0n_ || dataSize != tn_ * dipn_ * sdipn_) 
		sf_error ("%s: wrong dimensions modelSize=%d dataSize=%d",__FILE__, modelSize, dataSize);
	    break;
	case 2:
	    if (modelSize != tn_ * (xin_*sxin_ + dip0n_* sdip0n_) || dataSize != tn_ * dipn_ * sdipn_) 
		sf_error ("%s: wrong dimensions modelSize=%d dataSize=%d",__FILE__, modelSize, dataSize);
	    break;
    }

    sf_adjnull (adj, add, modelSize, dataSize, modl, data);

    pTableR = tableRefl_;
    pTableD = tableDiff_;

    for (sid = 0; sid < sdipn_; ++sid) { 
	    for (id = 0; id < dipn_; ++id) { 

			// diffraction part
			if (1 != invMod_) {	
				for (sixi = 0; sixi < sxin_; ++sixi) { 		
			    	for (ixi = 0; ixi < xin_; ++ixi) { 
						for (it = 0; it < tn_; ++it, ++pTableD) { 
						    curTime = to_ + it * td_;
						    t = *pTableD;		    
						    if (t > 0. && t < tLim_) {
								str [it] = t;
								tx[it] = 0.f;  // not sure if it is reasonable to calculate tx - may be too expensive
//							    tx [it] = anti * fabsf (p - p1) * dx;
								amp [it] = 1.;
						    } else {
								str[it] = curTime - 2.f * td_;
								tx[it]  = 0.f;
								amp[it] = 0.f;
						    }
						}	

						if (isAA_) {
						    sf_aastretch_define (str, tx, amp);
						    sf_chain (sf_halfint_lop, sf_aastretch_lop,
								      adj, true, tn_, tn_, tn_, modl + (sixi * xin_ + ixi) * tn_, data + (sid * dipn_ + id) * tn_, tmp);
						} else {
						    stretch_define (str);
						    stretch_lop (adj, true, tn_, tn_, modl + (sixi * xin_ + ixi) * tn_, data + (sid * dipn_ + id) * +tn_);
						}
		    		}
				}		
			}
	
			// reflection part
			if (0 != invMod_) {
			    offset = (1 == invMod_) ? 0 : xin_*sxin_*tn_;
	
				for (sid0 = 0; sid0 < sdip0n_; ++sid0) { 
				    for (id0 = 0; id0 < dip0n_; ++id0) { 
						for (it = 0; it < tn_; ++it, ++pTableR) {		
						    curTime = to_ + it * td_;
						    t = *pTableR;
						    if (t > 0. && t < tLim_) {
								str[it] = t;
								tx[it] = 0.f;  // not sure if it is reasonable to calculate tx - may be too expensive
//								tx[it] = fabsf(anti*sx)/t;
								amp[it] = 1.f;
						    } else {
								str[it] = curTime - 2.f * td_;
								tx[it]  = 0.f;
								amp[it] = 0.f;
						    }		
						}
					
						if (isAA_) {
						    sf_aastretch_define (str, tx, amp);
						    sf_chain (sf_halfint_lop, sf_aastretch_lop,
								      adj, true, tn_, tn_, tn_, modl + (sid0*dip0n_ + id0) * tn_ + offset, data + (sid*dipn_+id) * tn_, tmp);
						} else {
						    stretch_define (str);
						    stretch_lop (adj, true, tn_, tn_, modl + (sid0*dip0n_ + id0) * tn_ + offset, data + (sid*dipn_+id) * tn_);
						}
			    	}
				}	
			}
    	}
	}	
	
   	return;
}
