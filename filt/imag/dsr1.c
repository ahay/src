/* 2-D prestack migration/modeling by split-step DSR */
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

#include "dsr1.h"

void dsr1 (bool verb                  /* verbosity flag */, 
	   bool inv                   /* migration/modeling flag */, 
	   float eps                  /* stability factor */,  
	   int nw, float dw           /* frequency */, 
	   int nz, float dz           /* depth */,
	   int nx, float dx           /* midpoint */,
	   int nh, float dh, float h0 /* offset */,
	   float **vt                 /* slowness [nx][nz] */, 
	   float *v                   /* reference slowness [nz] */,
	   float complex ***cp        /* data [nx][nh][nw] */, 
	   float **q                  /* image [nx][nz] */)
/*< Apply migration/modeling >*/
{
    int nk,ng,iz,iw,ix,jx,ih,jh;
    float k,dk,fk, g,dg,fg, s,r;
    float complex cshift,**pp, w2;
    kiss_fft_cfg xforw, xinvs, hforw, hinvs;
    
    /* determine wavenumber sampling, pad by 2 */
    nk = nx*2;
    nk *= 2;
    dk = 2.0*SF_PI/(nk*dx);
    fk = -SF_PI/dx;
    
    ng = nh*2;
    nh *= 2;
    dg = 2.0*SF_PI/(ng*dh);
    fg = -SF_PI/dh;

    xforw = kiss_fft_alloc(nk,0,NULL,NULL);
    xinvs = kiss_fft_alloc(nk,1,NULL,NULL);
    hforw = kiss_fft_alloc(nk,0,NULL,NULL);
    hinvs = kiss_fft_alloc(nk,1,NULL,NULL);

    if (NULL == xforw || NULL == xinvs || 
	NULL == hforw || NULL == xinvs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    /* allocate workspace */
    pp  = sf_complexalloc2 (ng,nk);
  
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
		for (ih=0; ih < nh; ih++) {
		    pp[ix][ih] = q[ix][nz-1];
		}
		for (ih=nh; ih < ng; ih++) {
		    pp[ix][ih] = 0.;
		}
	    }
	    for (ix=nx; ix<nk; ix++) {
		for (ih=nh; ih < ng; ih++) {
		    pp[ix][ih] = 0.;
		}
	    }

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		for (ih=0; ih<ng; ih++) {
		    kiss_fft_stride(xforw,(kiss_fft_cpx *) (pp[0]+ih),
				    (kiss_fft_cpx *) (pp[0]+ih),ng);
		}		
		for (ix=0; ix<nk; ix++) {
		    kiss_fft(hforw,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);

		    jx = (ix < nk/2)? ix + nk/2: ix - nk/2;
		    k = fk + jx*dk;
		    for (ih=0; ih<ng; ih++) {
			jh = (ih < ng/2)? ih + ng/2: ih - ng/2;
			g = fg + jh*dg;
	
			s = 0.5*(k-g);
			r = 0.5*(k+g);
			s *= s;
			r *= r;

			cshift = cexpf((2.*csqrtf(w2*v[iz])-
					csqrtf(w2*v[iz]+s) -
					csqrtf(w2*v[iz]+r))*dz); 
			pp[ix][ih] *= cshift/(nk*ng); 
			/* FFT scaling included */
		    }
 	
		    kiss_fft(hinvs,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}
		for (ih=0; ih<ng; ih++) {
		    kiss_fft_stride(xinvs,(kiss_fft_cpx *) (pp[0]+ih),
				    (kiss_fft_cpx *) (pp[0]+ih),ng);
		}

		for (ix=0; ix<nx; ix++) {
		    cshift = cexpf(-csqrtf(w2*vt[ix][iz])*dz);
		    for (ih=0; ih < nh; ih++) {
			pp[ix][ih] = q[ix][iz] + pp[ix][ih]*cshift;
		    }
		    for (ih=nh; ih < ng; ih++) {
			pp[ix][ih] = 0.;
		    }
		}
		for (ix=nx; ix<nk; ix++) {
		    for (ih=0; ih < ng; ih++) {
			pp[ix][ih] = 0.;
		    }
		}
	    }

	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    cp[ix][ih][iw] = pp[ix][ih];
		}
	    }
	} else { /* migration */
    
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    pp[ix][ih] = cp[ix][ih][iw];
		}
	    }

	    /* loop over migrated depths z */
	    for (iz=0; iz<nz-1; iz++) {
		/* accumulate image (summed over frequency) */
		for (ix=0; ix<nx; ix++) { 
		    for (ih=0; ih < nh; ih++) {
			q[ix][iz] += crealf(pp[ix][ih]);
		    }
		}
		
		for (ix=0; ix<nx; ix++) {
		    cshift = conjf(cexpf(-csqrtf(w2*vt[ix][iz])*dz));
		    for (ih=0; ih < nh; ih++) {
			pp[ix][ih] *= cshift;
		    }
		    for (ih=nh; ih < ng; ih++) {
			pp[ix][ih] = 0.;
		    }
		}

		for (ix=nx; ix<nk; ix++) {
		    for (ih=0; ih < ng; ih++) {
			pp[ix][ih] = 0.;
		    }
		}

		for (ih=0; ih<ng; ih++) {
		    kiss_fft_stride(xforw,(kiss_fft_cpx *) (pp[0]+ih),
				    (kiss_fft_cpx *) (pp[0]+ih),ng);
		}

		for (ix=0; ix<nk; ix++) {
		    kiss_fft(hforw,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);

		    jx = (ix < nk/2)? ix + nk/2: ix - nk/2;
		    k = fk + jx*dk;
		    for (ih=0; ih<ng; ih++) {
			jh = (ih < ng/2)? ih + ng/2: ih - ng/2;
			g = fg + jh*dg;
	
			s = 0.5*(k-g);
			r = 0.5*(k+g);
			s *= s;
			r *= r;
	  
			cshift = conjf(cexpf((2.*csqrtf(w2*v[iz])-
					      csqrtf(w2*v[iz]+s) -
					      csqrtf(w2*v[iz]+r))*dz)); 

			pp[ix][ih] *= cshift/(nk*ng);
			/* Fourier scaling included */
		    }

		    kiss_fft(hinvs,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}
		
		for (ih=0; ih<ng; ih++) {
		    kiss_fft_stride(xinvs,(kiss_fft_cpx *) (pp[0]+ih),
				    (kiss_fft_cpx *) (pp[0]+ih),ng);
		}
	    }

	    for (ix=0; ix<nx; ix++) { 
		for (ih=0; ih < nh; ih++) {
		    q[ix][nz-1] += crealf(pp[ix][ih]);
		}
	    }
	}
    }

    free (*pp);
    free (pp);
    free (xforw);
    free (xinvs);
    free (hforw);
    free (hinvs);
}

/* 	$Id: split1.c 790 2004-09-10 18:51:51Z fomels $	 */
