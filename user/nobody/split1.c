/* 2-D split-step migration/modeling */
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

#include "split1.h"

void split1 (bool verb          /* verbosity flag */, 
	     bool inv           /* migration/modeling flag */, 
	     float eps          /* stability factor */,  
	     int nw, float dw   /* frequency */, 
	     int nz, float dz   /* depth */,
	     int nx, float dx   /* midpoint */,
	     float **vt         /* slowness [nx][nz] */, 
	     float *v           /* reference slowness [nz] */,
	     float complex **cp /* data [nx][nw] */, 
	     float **q          /* image [nx][nz] */)
/*< Apply migration/modeling >*/
{
    int nk,iz,iw,ix,jx;
    float k,dk,fk;
    float complex cshift,*pp,w2;
    kiss_fft_cfg forw, invs;

    /* determine wavenumber sampling, pad by 2 */
    nk = nx*2;
    nk *= 2;
    dk = 2.0*SF_PI/(nk*dx);
    fk = -SF_PI/dx;

    forw = kiss_fft_alloc(nk,0,NULL,NULL);
    invs = kiss_fft_alloc(nk,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    /* allocate workspace */
    pp  = sf_complexalloc (nk);
  
    if (!inv) { /* prepare image for migration */
	for (ix=0; ix<nx; ix++) {      
	    for (iz=0; iz<nz; iz++)
		q[ix][iz] = 0.0;
	}
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w2 = dw*(eps+iw*I);
	w2 *= w2;

	if (inv) { /* modeling */
	    for (ix=0; ix<nx; ix++) {
		pp[ix] = q[ix][nz-1];
	    }
	    for (ix=nx; ix<nk; ix++) {
		pp[ix] = 0.;
	    }

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		kiss_fft(forw,(const kiss_fft_cpx *) pp, 
			 (kiss_fft_cpx *) pp);

		for (ix=0; ix<nk; ix++) {
		    jx = (ix < nk/2)? ix + nk/2: ix - nk/2;
		    k = fk + jx*dk;
		    k *= k;
	  
		    cshift = cexpf((csqrtf(w2*v[iz])-
				    csqrtf(w2*v[iz]+k))*dz); 
		    pp[ix] *= cshift/nk; 
		    /* FFT scaling included */
		}
		
		kiss_fft(invs,(const kiss_fft_cpx *) pp, 
			 (kiss_fft_cpx *) pp);

		for (ix=0; ix<nx; ix++) {
		    cshift = cexpf(-csqrtf(w2*vt[ix][iz])*dz);
		    pp[ix]= q[ix][iz] + pp[ix]*cshift;
		}
		for (ix=nx; ix<nk; ix++) {
		    pp[ix] = 0.;
		}
	    }

	    for (ix=0; ix<nx; ix++) {
		cp[ix][iw] = pp[ix];
	    }
	} else { /* migration */
    
	    for (ix=0; ix<nx; ix++) {
		pp[ix] = cp[ix][iw];
	    }
	    for (ix=nx; ix<nk; ix++) {
		pp[ix] = 0.;
	    }

	    /* loop over migrated depths z */
	    for (iz=0; iz<nz-1; iz++) {
		/* accumulate image (summed over frequency) */
		for (ix=0; ix<nx; ix++) { 
		    q[ix][iz] += crealf(pp[ix]);
		}
		
		for (ix=0; ix<nx; ix++) {
		    cshift = conjf(cexpf(-0.5*csqrtf(w2*vt[ix][iz])*dz));
		    pp[ix] *= cshift;
		}
		for (ix=nx; ix<nk; ix++) {
		    pp[ix] = 0.;
		}

		kiss_fft(forw,(const kiss_fft_cpx *) pp, 
			 (kiss_fft_cpx *) pp);

		for (ix=0; ix<nk; ix++) {
		    jx = (ix < nk/2)? ix + nk/2: ix - nk/2;
		    k = fk + jx*dk;
		    k *= k;
	  
		    cshift = conjf(cexpf((0.5*(csqrtf(w2*v[iz])+
					       csqrtf(w2*v[iz+1]))-
					  csqrtf(w2*v[iz]+k))*dz));
		    pp[ix] *= cshift/nk;
		    /* Fourier scaling included */
		}

		kiss_fft(invs,(const kiss_fft_cpx *) pp, 
			 (kiss_fft_cpx *) pp);
	
		for (ix=0; ix<nx; ix++) {
		    cshift = conjf(cexpf(-0.5*csqrtf(w2*vt[ix][iz+1])*dz));
		    pp[ix] *= cshift;
		}
	    }

	    for (ix=0; ix<nx; ix++) { 
		q[ix][nz-1] += crealf(pp[ix]);
	    }
	}
    }
    
    free (pp);
    free (forw);
    free (invs);
}

/* 	$Id$	 */
