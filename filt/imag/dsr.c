/* DSR modeling/migration in v(z) */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>

#include <rsf.h>

#include "dsr.h"
#include "pshift.h"

static float eps, *vt, dw, dz, da;
static int nz, nw, na;
static float complex *pp;
static kiss_fftr_cfg forw, invs;

void dsr_init (float eps1 /* regularization */, 
	       int nt     /* time samples */, 
	       float dt   /* time sampling */, 
	       int nz1    /* depth samples */, 
	       float dz1  /* depth sampling */, 
	       float *vt1 /* velocity/slowness */, 
	       bool depth /* depth or time migration */,
	       char rule   /* interpolation rule */,
	       int na1     /* angle samples */,
	       float da1   /* angle sampling */)
/*< initialize >*/
{
    eps = eps1;     
    nz = nz1; 
    dz = dz1;
    vt = vt1;
    na = na1;
    da = da1;

    /* determine frequency sampling */
    nw = nt/2+1;
    dw = 2.0*SF_PI/(nt*dt);
    
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);
    pshift_init(depth,rule);
} 

void dsr_close ()
/*< free workspace >*/
{    
    free (pp);  
    free (forw);
    free (invs);
}

void dsr (bool inv /* modeling or migration */, 
	  float kx /* midpoint wavenumber */, 
	  float kh /* half-offset wavenumber */, 
	  float *p /* time trace */, 
	  float **q /* depth trace */)
/*< apply >*/
{
    int iz,iw;
    float s, r, a;
    float complex w, k;

    s = 0.5*(kx-kh);
    r = 0.5*(kx+kh);
    s *= s;
    r *= r;

    if (inv) { /* modeling */
	for (iw=0; iw<nw; iw++) {
	    pp[iw] = q[nz-1][0];
	}
	
	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w = (eps + I*iw)*dw;
		pp[iw] = q[iz][0] + pp[iw]*
		    cexpf(-(pshift(w,r,vt[iz],vt[iz+1])+
			    pshift(w,s,vt[iz],vt[iz+1]))*dz);
	    }
	}

	kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    } else { /* migration */
	kiss_fftr(forw, p, (kiss_fft_cpx *) pp);

	/* loop over migrated times z */
	for (iz=0; iz<nz-1; iz++) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w = (eps+I*iw)*dw;

		/* find angle */
		
		/* accumulate image (summed over frequency) */
		q[iz][0] += crealf(pp[iw]);
		pp[iw] *= conjf(cexpf(-(pshift(w,r,vt[iz],vt[iz+1])+
					pshift(w,s,vt[iz],vt[iz+1]))*dz));
	    }
	}

	for (iw=0; iw<nw; iw++) {
	    q[nz-1][0] += crealf(pp[iw]);
	}
    }
}

/* 	$Id$	 */

