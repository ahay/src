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
#include "pshift.h"

static float eps, *vt, *vz, *eta, dw, dz;
static int nz, nw;
static sf_complex *pp;
static kiss_fftr_cfg forw, invs;

void gazdag_init (float eps1  /* regularization */, 
		  int nt      /* time samples */, 
		  float dt    /* time sampling */, 
                  int nz1     /* depth samples */, 
		  float dz1   /* depth sampling */, 
		  float *vt1  /* velocity (time) or slowness (depth) */, 
		  float *vz1,
		  float *eta1,
		  bool depth  /* depth (or time) */,
		  char rule   /* interpolation rule */)
/*< Initialize >*/
{
    eps = eps1; 
    nz = nz1; 
    dz = dz1;
    vt = vt1;
    vz = vz1;
    eta = eta1;

    /* determine frequency sampling */
    nw = nt/2+1;
    dw = 2.0*SF_PI/(nt*dt);
    
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);
    pshift_init(depth,rule);
}

void gazdag_close ()
/*< Free allocated storage >*/
{    
    free (pp);  
    free (forw);
    free (invs);
}

void gazdag (bool inv  /* modeling (or migration) */,
	     float k2  /* wavenumber squared */, 
	     float *p  /* data [nt] */, 
	     float *q  /* image [nz] */)
/*< Run >*/
{
    int iz,iw;
    sf_complex w2, kz;
        
    if (inv) { /* modeling */
        for (iw=0; iw<nw; iw++) {
            pp[iw] = sf_cmplx(q[nz-1],0.);
        }

        /* loop over migrated times z */
        for (iz=nz-2; iz>=0; iz--) {
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
                w2 = sf_cmplx(eps*dw,iw*dw);
		kz = pshift(w2,k2,vt[iz],vt[iz+1],vz[iz],eta[iz]);
#ifdef SF_HAS_COMPLEX_H
                pp[iw] = pp[iw]*cexpf(-kz*dz) + q[iz];
#else
		pp[iw] = sf_cadd(sf_cmul(pp[iw],cexpf(sf_crmul(kz,-dz))),
				 sf_cmplx(q[iz],0.));
#endif
            }
        }

        kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    } else { /* migration */
	kiss_fftr(forw, p, (kiss_fft_cpx *) pp);
    
        /* loop over migrated times z */
        for (iz=0; iz<nz-1; iz++) {
            /* initialize migrated sample */
            q[iz] = 0.0;
      
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
		w2 = sf_cmplx(eps*dw,iw*dw);

                /* accumulate image (summed over frequency) */
                q[iz] += crealf(pp[iw]);
		kz = pshift(w2,k2,vt[iz],vt[iz+1],vz[iz],eta[iz]);
#ifdef SF_HAS_COMPLEX_H
                pp[iw] *= conjf(cexpf(-kz*dz));
#else
		pp[iw] = sf_cmul(pp[iw],conjf(cexpf(sf_crmul(kz,-dz))));
#endif
            }
        }

	q[nz-1] = 0.;
	for (iw=0; iw<nw; iw++) {
	    q[nz-1] += crealf(pp[iw]);
	}
    }
}

/*      $Id: gazdag.c 7107 2011-04-10 02:04:14Z ivlad $     */
