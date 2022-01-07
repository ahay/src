/* Digital oclet transform in f-k domain */
/*
  Copyright (C) 2009 University of Texas at Austin
   
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

#include "fkoclet3.h"

static int h, nha;
static float dh, da, kx, ky, h0, a0, w, epsilon, maxe;
static bool inv, unit, dwt, amplitude;
static sf_complex *t, t1, t2;
static float *wei;
static void (*transform)(bool);

static sf_complex fkocpredict(bool forw, sf_complex tt, int i, int j)
/* Predict forward or backward */
{
    if (dwt) {
	return(tt);
    } else {
	float h1, h2, eps1, amp1, eps2, amp2, amp, phase1, phase2, k1, k2;
	sf_complex oper;
	
	if (forw) {
	    h1 = h0 + (i)*dh;
	    h2 = h0 + (i+j)*dh;

	    k1 = kx*cosf(a0+(i)*da)+ky*sinf(a0+(i)*da);
	    k2 = kx*cosf(a0+(i+j)*da)+ky*sinf(a0+(i+j)*da);
	} else {
	    h1 = h0 + (i+j)*dh;
	    h2 = h0 + (i)*dh;

	    k1 = kx*cosf(a0+(i+j)*da)+ky*sinf(a0+(i+j)*da);
	    k2 = kx*cosf(a0+(i)*da)+ky*sinf(a0+(i)*da);
	}

	if (fabsf (w) > FLT_EPSILON) {
	    eps1 = 2.*fabsf(k1*h1/w);
	    eps2 = 2.*fabsf(k2*h2/w);
	    if (eps1 <= maxe && eps2 <= maxe) {
		eps1 = hypotf (1.,eps1);
		eps2 = hypotf (1.,eps2);
                 
		amp1 = 1./eps1+eps1;
		amp2 = 1./eps2+eps2;
		phase1 = 1-eps1+logf(0.5*(1.+eps1));
		phase2 = 1-eps2+logf(0.5*(1.+eps2));
		
		if (amplitude) {
		    amp = expf(0.5*(eps1-logf(amp1)+logf(amp2)-eps2)); 
		} else {
		    amp = 1.;
		}
		
		phase1 = SF_PI*(phase2-phase1)*w;
		oper = sf_cmplx(amp*cosf(phase1),amp*sinf(phase1));
#ifdef SF_HAS_COMPLEX_H
		tt = tt*oper;
#else
		tt = sf_cmul(tt,oper);
#endif
	    } 
		//else {
		//tt = tt;
	    //}
	} 
	//else {
	//    tt = tt;
	//}
	return(tt);
    }
}

static void fkochaar(bool adj) 
/* Lifting Haar transform in offset */
{
    int i, j;

    if (adj) {
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		if (inv) {
		    t2 = t[i+j];
		    t2 = fkocpredict(false,t2,i,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i] -= t2/2;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(t2,-0.5));
#endif
		    t1 = t[i];
		    t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += t1;
#else
		    t[i+j] = sf_cadd(t[i+j],t1);		    
#endif
		} else {
		    t1 = t[i];
		    t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += t1/2;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(t1,0.5));
#endif
		    t2 = t[i+j];
		    t2 = fkocpredict(false,t2,i,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i] -= t2;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(t2,-1.));
#endif
		}
	    }
	}
    } else {
	for (j=1; j <= h/2; j *= 2) {
	    for (i=0; i < h-j; i += 2*j) {


		t1 = t[i];
		t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H
		t[i+j] -= t1;
#else
		t[i+j] = sf_cadd(t[i+j],sf_crmul(t1,-1.));
#endif
		t2 = t[i+j];
		t2 = fkocpredict(false,t2,i,j);
#ifdef SF_HAS_COMPLEX_H
		t[i] += t2/2;
#else
		t[i] = sf_cadd(t[i],sf_crmul(t2,0.5));
#endif
	    }	    
	}
    }
}

static void fkoclinear(bool adj) 
/* Lifting linear-interpolation transform in offset */
{
    int i, j;

    if (adj) {
	for (j=h/2; j >= 1; j /= 2) {
	    if (inv) {
		for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i-j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i]   -= (t1+t2)/4;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),-0.25));
#endif
		}
		t1 = t[j];
		t1 = fkocpredict(false,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
		t[0] -= t1/2;
#else
		t[0] = sf_cadd(t[0],sf_crmul(t1,-0.5));
#endif
		for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i+2*j];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t1+t2)/2;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),0.5));
#endif
		}
		t1 = t[i];
		t1 = fkocpredict(true,t1,i,j);	
#ifdef SF_HAS_COMPLEX_H	
		if (i+j < h) t[i+j] += t1;
#else
		if (i+j < h) t[i+j] = sf_cadd(t[i+j],t1);
#endif
	    } else {
		for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H	
		    t[i+j] += t1/4;
		    t[i-j] += t2/4;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(t1,0.25));
		    t[i-j] = sf_cadd(t[i-j],sf_crmul(t2,0.25));
#endif
		}
		t1 = t[0];
		t1 = fkocpredict(true,t1,0,j);
#ifdef SF_HAS_COMPLEX_H	
		t[j] += t1/2;
#else
		t[j] = sf_cadd(t[j],sf_crmul(t1,0.5));
#endif
		for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i+j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H	
		    t[i]     -= t1/2;
		    t[i+2*j] -= t2/2;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(t1,-0.5));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(t2,-0.5));
#endif
		}
		t1 = t[i+j];
		t1 = fkocpredict(false,t1,i,j);	 
#ifdef SF_HAS_COMPLEX_H	
		if (i+j < h) t[i] -= t1;
#else
		if (i+j < h) t[i] = sf_cadd(t[i],sf_cneg(t1));
#endif
	    }
	}
    } else {
	for (j=1; j <= h/2; j *= 2) {
	    for (i=0; i < h-2*j; i += 2*j) {
		t1 = t[i];
		t2 = t[i+2*j];
		t1 = fkocpredict(true,t1,i,j);
		t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H	
		t[i+j] -= (t1+t2)/2;
		/* d = o - P e */
#else
		t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),-0.5));
#endif
	    }
	    t1 = t[i];
	    t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H
	    if (i+j < h) t[i+j] -= t1;    
#else
	    if (i+j < h) t[i+j] = sf_cadd(t[i+j],sf_cneg(t1));
#endif

	    t2 = t[j];
	    t2 = fkocpredict(false,t2,0,j);	 
#ifdef SF_HAS_COMPLEX_H
	    t[0] += t2/2;
#else
	    t[0] = sf_cadd(t[0],sf_crmul(t2,0.5));
#endif
	    for (i=2*j; i < h-j; i += 2*j) {
		t1 = t[i+j];
		t2 = t[i-j];
		t1 = fkocpredict(false,t1,i,j);
		t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t1+t2)/4;
		/* s = e + U d */
#else
		t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),0.25));
#endif
	    }
	}
    }
}

static void fkocbior(bool adj) 
/* Lifting CDF 9/7 biorthogonal wavelet transform in offset */
{
    int i, j;
    float a;

    if (adj) {
	for (j=h/2; j >= 1; j /= 2) {
	    if (inv) {                            /*reverse dwt9/7 transform*/
                a= 1.230174105f;
	        for (i=2*j; i < h-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]  /= a;
#else
	            t[i] = sf_crmul(t[i],(1/a));
#endif
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] /= a;                         /*left boundary*/
#else
		t[0] = sf_crmul(t[0],(1/a));
#endif
	        for (i=0; i < h-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] *= a;
#else
	            t[i+j] = sf_crmul(t[i+j],a);
#endif
		    /* Undo Scale */
	        }
#ifdef SF_HAS_COMPLEX_H
	        if (i+j < h) t[i+j] *= a;         /*right boundary*/  
#else
	        if (i+j < h) t[i+j] = sf_crmul(t[i+j],a);
#endif
                a= -0.4435068522f;
	        for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i-j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i]   += (t1+t2)*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),a));
#endif
		    /* Undo Update 2 */
	        }
		t1 = t[j];
		t1 = fkocpredict(false,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	        t[0] += 2*a*t1;                  /*left boundary*/
#else
		t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(t1,a),2));
#endif
	        a = -0.8829110762f;
	        for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i+2*j];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t1+t2)*a;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),a));
#endif
		    /* Undo Predict 2 */
	        }
		t1 = t[i];
		t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < h) t[i+j] += 2*a*t1;  /*right boundary*/  
#else
		if (i+j < h) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(t1,a),2));
#endif
                    /* Undo Step 2 */

                a= 0.05298011854f;
	        for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i-j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i]   += (t1+t2)*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),a));
#endif
		    /* Undo Update 1 */
	        }
		t1 = t[j];
		t1 = fkocpredict(false,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	        t[0] += 2*a*t1;                  /*left boundary*/
#else
		t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(t1,a),2));
#endif
	        a = 1.586134342f;
	        for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i+2*j];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t1+t2)*a;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),a));
#endif
		    /* Undo Predict 1 */
	        }
		t1 = t[i];
		t1 = fkocpredict(true,t1,i,j);	 
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < h) t[i+j] += 2*a*t1;  /*right boundary*/  
#else
		if (i+j < h) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(t1,a),2));
#endif
                    /* Undo Step 1 */
	    } else {                               /*adjoint transform*/
                a= 1.230174105f;
	        for (i=2*j; i < h-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]  *= a;
#else
	            t[i] = sf_crmul(t[i],a);
#endif
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] *= a;                         /*left boundary*/
#else
		t[0] = sf_crmul(t[0],a);
#endif
	        for (i=0; i < h-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] /= a;
#else
	            t[i+j] = sf_crmul(t[i+j],(1/a));
#endif
		    /* Undo Scale */
	        }
#ifdef SF_HAS_COMPLEX_H
	        if (i+j < h) t[i+j] /= a;         /*right boundary*/  
#else
	        if (i+j < h) t[i+j] = sf_crmul(t[i+j],(1/a));
#endif
                a= -0.4435068522f;
	        for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j]   -= t1*a;
                    t[i-j]   -= t2*a;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(t1,(-a)));
		    t[i-j] = sf_cadd(t[i-j],sf_crmul(t2,(-a)));
#endif
		    /* Undo Update 2 */
	        }
		t1 = t[0];
		t1 = fkocpredict(true,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	        t[j] -= 2*a*t1;                  /*left boundary*/
#else
		t[j] = sf_cadd(t[j],sf_crmul(sf_crmul(t1,a),-2));
#endif
	        a = -0.8829110762f;
	        for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i+j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
                    t[i]     -= t1*a;
                    t[i+2*j] -= t2*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(t1,(-a)));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(t2,(-a)));
#endif
		    /* Undo Predict 2 */
	        }
		t1 = t[i+j];
		t1 = fkocpredict(false,t1,i,j);
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < h) t[i] -= 2*a*t1;  /*right boundary*/  
#else
		if (i+j < h) t[i] = sf_cadd(t[i],sf_crmul(sf_crmul(t1,a),-2));
#endif
                    /* Undo Step 2 */
 
                a= 0.05298011854f;
	        for (i=2*j; i < h-j; i += 2*j) {
		    t1 = t[i];
		    t2 = t[i];
		    t1 = fkocpredict(true,t1,i,j);
		    t2 = fkocpredict(false,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		    t[i+j]   -= t1*a;
                    t[i-j]   -= t2*a;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(t1,(-a)));
		    t[i-j] = sf_cadd(t[i-j],sf_crmul(t2,(-a)));
#endif		    /* Undo Update 1 */
	        }
		t1 = t[0];
		t1 = fkocpredict(true,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	        t[j] -= 2*a*t1;                  /*left boundary*/
#else
		t[j] = sf_cadd(t[j],sf_crmul(sf_crmul(t1,a),-2));
#endif
	        a = 1.586134342f;
	        for (i=0; i < h-2*j; i += 2*j) {
		    t1 = t[i+j];
		    t2 = t[i+j];
		    t1 = fkocpredict(false,t1,i,j);
		    t2 = fkocpredict(true,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
                    t[i]     -= t1*a;
                    t[i+2*j] -= t2*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(t1,(-a)));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(t2,(-a)));
#endif
		    /* Undo Predict 1 */
	        }
		t1 = t[i+j];
		t1 = fkocpredict(false,t1,i,j);	 
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < h) t[i] -= 2*a*t1;  /*right boundary*/  
#else
		if (i+j < h) t[i] = sf_cadd(t[i],sf_crmul(sf_crmul(t1,a),-2));
#endif                    /* Undo Step 1 */
	    }
	}
    } else {
	for (j=1; j <= h/2; j *= 2) {        /*different scale*/
	    a = -1.586134342f;
	    for (i=0; i < h-2*j; i += 2*j) {
		t1 = t[i];
		t2 = t[i+2*j];
		t1 = fkocpredict(true,t1,i,j);
		t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
	  	 t[i+j] += (t1+t2)*a;
#else
	   	 t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),a));
#endif
		/* Predict 1 */
	    }	 
	    t1 = t[i];
	    t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H	 
            if (i+j < h) t[i+j] += 2*a*t1;  /*right boundary*/  
#else
	    if (i+j < h) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(t1,a),2));
#endif
            a= -0.05298011854f;
	    t1 = t[j];
	    t1 = fkocpredict(false,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	    t[0] += 2*a*t1;                  /*left boundary*/
#else
	    t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(t1,a),2));
#endif
	    for (i=2*j; i < h-j; i += 2*j) {
		t1 = t[i+j];
		t2 = t[i-j];
		t1 = fkocpredict(false,t1,i,j);
		t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t1+t2)*a;
#else
		t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),a));
#endif
		/* Update 1 */
	    }
                /* Step 1 */
	    a = 0.8829110762f;
	    for (i=0; i < h-2*j; i += 2*j) {
		t1 = t[i];
		t2 = t[i+2*j];
		t1 = fkocpredict(true,t1,i,j);
		t2 = fkocpredict(false,t2,i+j,j);
#ifdef SF_HAS_COMPLEX_H
	  	 t[i+j] += (t1+t2)*a;
#else
	   	 t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cadd(t1,t2),a));
#endif
		/* Predict 2 */
	    }
	    t1 = t[i];
	    t1 = fkocpredict(true,t1,i,j);
#ifdef SF_HAS_COMPLEX_H	 
            if (i+j < h) t[i+j] += 2*a*t1;  /*right boundary*/  
#else
	    if (i+j < h) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(t1,a),2));
#endif
            a= 0.4435068522f;
	    t1 = t[j];
	    t1 = fkocpredict(false,t1,0,j);
#ifdef SF_HAS_COMPLEX_H
	    t[0] += 2*a*t1;                  /*left boundary*/
#else
	    t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(t1,a),2));
#endif
	    for (i=2*j; i < h-j; i += 2*j) {
		t1 = t[i+j];
		t2 = t[i-j];
		t1 = fkocpredict(false,t1,i,j);
		t2 = fkocpredict(true,t2,i-j,j);
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t1+t2)*a;
#else
		t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(t1,t2),a));
#endif
		/* Update 2 */
	    }
                /* Step 2 */
            a= 1/(1.230174105f);
	    for (i=0; i < h-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i+j] *= a;
#else
	        t[i+j] = sf_crmul(t[i+j],a);
#endif
	    }	 
#ifdef SF_HAS_COMPLEX_H
	    if (i+j < h) t[i+j] *= a;         /*right boundary*/  
#else
	    if (i+j < h) t[i+j] = sf_crmul(t[i+j],a);
#endif
#ifdef SF_HAS_COMPLEX_H
	    t[0] /= a;                         /*left boundary*/
#else
	    t[0] = sf_crmul(t[0],(1/a));
#endif
	    for (i=2*j; i < h-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i]  /= a;
#else
	        t[i] = sf_crmul(t[i],(1/a));
#endif
		/* Scale */
	    }
	}
    }
}

void fkoclet_init(int nha_in /* data size */, 
		  bool inv1, 
		  bool unit1, 
		  bool dwt1,
		  bool amp1,
		  float e,
		  char type,
		  float maxe1) 
/*< allocate space >*/
{
    int i, j;
    float wi;

    inv = inv1;
    unit = unit1;
    dwt = dwt1;
    amplitude = amp1;
    nha = nha_in;
    epsilon = e;
    maxe = maxe1;

    for (h=1; h < nha; h *= 2) ;
    t = sf_complexalloc(h);
    
    switch(type) {
	case 'h': 
	    transform = fkochaar;
	    break;
	case 'l':
	    transform = fkoclinear;
	    break;
	case 'b':
	    transform = fkocbior;
	    break;
	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }

    if (unit) {
	wei = sf_floatalloc(h);

	wei[0] = sqrtf((float) h);
	wi = 0.5;	
	for (j=1; j <= h/2; j *= 2, wi *= 2) {
	    
	    for (i=0; i < h-j; i += 2*j) {
		wei[i+j] = sqrtf(wi);
	    }
	}
    }
}


void fkoclet_close(void) 
/*< deallocate space >*/
{
    free (t);
    if (unit) free(wei);
}

void fkoclet_lop(bool adj, bool add, int nx, int ny, 
		 sf_complex *x, sf_complex *y, float w_in, 
		 float kx_in, float ky_in, float h0_in, 
		 float dh_in, float a0_in, float da_in)
/*< linear operator >*/
{
    int it, i, j;
    w = w_in;
    kx = kx_in;
    ky = ky_in;
    h0 = h0_in;
    dh = dh_in;
    a0 = a0_in;
    da = da_in;

    sf_cadjnull (adj,add,nx,ny,x,y);

    if (adj) {
	for (it=0; it < ny; it++) {
	    t[it]=y[it];
	}
	for (it=ny; it < h; it++) {
	    t[it] = sf_cmplx(0.,0.);
	}
    } else{
	t[0] = x[0];
	it = 1;
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		if (it < nx) {
		    t[i+j]=x[it];
		    it++;
		} else {
		    t[i+j]=sf_cmplx(0.,0.);
		}
	    }	    	    
	}

	if (unit) {
	    for (it=0; it < h; it++) {
		if (inv) {
#ifdef SF_HAS_COMPLEX_H
		    t[it] /= wei[it];
#else
		    t[it] = sf_crmul(t[it],1.0f/wei[it]);
#endif
		} else {
#ifdef SF_HAS_COMPLEX_H
		    t[it] *= wei[it];
#else
		    t[it] = sf_crmul(t[it],wei[it]);
#endif
		}
	    }
	}
    } 

    transform((bool)!adj);    

    if (adj) {
	if (unit) {
	    for (it=0; it < h; it++) {
#ifdef SF_HAS_COMPLEX_H
		t[it] *= wei[it];
#else
		t[it] = sf_crmul(t[it],wei[it]);
#endif
	    }
	}

#ifdef SF_HAS_COMPLEX_H
	x[0] += t[0];
#else
	x[0] = sf_cadd(x[0],t[0]);
#endif
	it = 1;
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		x[it] += t[i+j];
#else
		x[it] = sf_cadd(x[it],t[i+j]);
#endif
		it++;
		if (it >= nx) return;
	    }	    	    
	}
    } else {
	for (it=0; it < ny; it++) {
#ifdef SF_HAS_COMPLEX_H
	    y[it] += t[it];
#else
	    y[it] = sf_cadd(y[it],t[it]);
#endif
	}
    } 
}
/* 	$Id$	 */


