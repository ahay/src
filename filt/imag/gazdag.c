/* Gazdag zero-offset phase-shift migration/modeling */
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
/*^*/

#include "gazdag.h"

static float eps, dz, *vt, dw;
static int nz, nw;
static float complex *pp;
static bool depth;
static kiss_fftr_cfg forw, invs;

void gazdag_init (float eps1  /* regularization */, 
		  int nt      /* time samples */, 
		  float dt    /* time sampling */, 
                  int nz1     /* depth samples */, 
		  float dz1   /* depth sampling */, 
		  float *vt1  /* velocity (time) or slowness (depth) */, 
		  bool depth1 /* depth (or time) */)
/*< Initialize >*/
{
    eps = eps1; 
    nz = nz1; dz = dz1;
    vt = vt1;
    depth = depth1;

    /* determine frequency sampling */
    nw = nt/2+1;
    dw = 2.0*SF_PI/(nt*dt);
    
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);
}

void gazdag_close ()
/*< Free allocated storage >*/
{    
    free (pp);  
    free (forw);
    free (invs);
}

void gazdag (bool inv /* modeling (or migration) */, 
	     float k2 /* wavenumber squared */, 
	     float *p /* data [nt] */, 
	     float *q /* image [nz] */)
/*< Run >*/
{
    int iz,iw;
    float complex cshift, w2;
        
    if (inv) { /* modeling */
        for (iw=0; iw<nw; iw++) {
            pp[iw] = q[nz-1];
        }

        /* loop over migrated times z */
        for (iz=nz-2; iz>=0; iz--) {
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
                w2 = (eps + I*iw)*dw;

                if (depth) {
                    w2 = w2*w2 * vt[iz] + k2;
                } else {
                    w2 = w2*w2 + vt[iz] * k2;
                }
        
                cshift = cexpf(-csqrtf(w2)*dz);
                pp[iw] = pp[iw]*cshift + q[iz];
            }
        }

        kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    } else { /* migration */
	kiss_fftr(forw, p, (kiss_fft_cpx *) pp);
    
        /* loop over migrated times z */
        for (iz=0; iz<nz; iz++) {
            /* initialize migrated sample */
            q[iz] = 0.0;
      
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
                /* accumulate image (summed over frequency) */
                q[iz] += crealf(pp[iw]);

                w2 = (eps + I*iw)*dw;

                if (depth) { /* depth migration */
                    w2 = w2*w2 * vt[iz] + k2;
                } else { /* time migration */
                    w2 = w2*w2 + vt[iz] * k2;
                }
        
                /* extrapolate down one migrated time step */
                cshift = cexpf(-csqrtf(w2)*dz);
                pp[iw] *= conjf(cshift);
            }
        }
    }
}

/*      $Id$     */
