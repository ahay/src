/* 3-D split-step migration/modeling */
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

#include "split3.h"

void split3 (bool verb           /* verbosity flag */, 
	     bool inv            /* migration/modeling flag */, 
	     float eps           /* stability factor */,  
	     int nw, float dw    /* frequency */, 
	     int nz, float dz    /* depth */,
	     int nx, float dx    /* cross-line midpoint */,
	     int ny, float dy    /* in-line midpoint */,
	     float ***vt         /* slowness [nx][ny][nz] */, 
	     float *v            /* reference slowness [nz] */,
	     float complex ***cp /* data [nx][ny][nw] */, 
	     float ***q          /* image [nx][ny][nz] */)
/*< Apply migration/modeling >*/
{
    int nkx,nky,iz,iw,ix,jx,iy,jy;
    float kx,dkx,fkx,ky,dky,fky;
    float complex cshift,**pp,w2, *ctrace;
    kiss_fft_cfg xforw, xinvs, yforw, yinvs;

    /* determine wavenumber sampling, pad by 2 */
    nkx = nx*2;
    nkx *= 2;
    dkx = 2.0*SF_PI/(nkx*dx);
    fkx = -SF_PI/dx;

    nky = ny*2;
    nky *= 2;
    dky = 2.0*SF_PI/(nky*dy);
    fky = -SF_PI/dy;

    xforw = kiss_fft_alloc(nkx,0,NULL,NULL);
    xinvs = kiss_fft_alloc(nkx,1,NULL,NULL);
    yforw = kiss_fft_alloc(nky,0,NULL,NULL);
    yinvs = kiss_fft_alloc(nky,1,NULL,NULL);

    if (NULL == xforw || NULL == xinvs || NULL == yforw || NULL == yinvs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    /* allocate workspace */
    pp  = sf_complexalloc2 (nky,nkx);
    ctrace = sf_complexalloc(nkx);
  
    if (!inv) { /* prepare image for migration */
	for (ix=0; ix<nx; ix++) {      
	    for (iy=0; iy<ny; iy++) {
		for (iz=0; iz<nz; iz++) {
		    q[ix][iy][iz] = 0.0;
		}
	    }
	}
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w2 = dw*(eps+iw*I);
	w2 *= w2;

	if (inv) { /* modeling */
	    for (ix=0; ix<nx; ix++) {
		for (iy=0; iy<ny; iy++) {
		    pp[ix][iy] = q[ix][iy][nz-1];
		}
		for (iy=ny; iy<nky; iy++) {
		    pp[ix][iy] = 0.;
		}
	    }
	    for (ix=nx; ix<nkx; ix++) {
		for (iy=0; iy<nky; iy++) {
		    pp[ix][iy] = 0.;
		}
	    }
	    
	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		for (iy=0; iy < nky; iy++) {
		    kiss_fft_stride(xforw,(const kiss_fft_cpx *) (pp[0]+iy), 
				    (kiss_fft_cpx *) ctrace,nky);
		    for (ix=0; ix<nkx; ix++) {
			pp[ix][iy] = ctrace[ix];
		    }
		}
		for (ix=0; ix < nkx; ix++) {
		    kiss_fft(yforw,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}

		for (ix=0; ix<nkx; ix++) {
		    jx = (ix < nkx/2)? ix + nkx/2: ix - nkx/2;
		    kx = fkx + jx*dkx;
		    kx *= kx;
		    
		    for (iy=0; iy<nky; iy++) {
			jy = (iy < nky/2)? iy + nky/2: iy - nky/2;
			ky = fky + jy*dky;
			ky = ky*ky + kx;
	  
			cshift = cexpf((csqrtf(w2*v[iz])-
					csqrtf(w2*v[iz]+ky))*dz); 
			pp[ix][iy] *= cshift/(nkx*nky); 
			/* FFT scaling included */
		    }
		}
		    
		for (ix=0; ix < nkx; ix++) {
		    kiss_fft(yinvs,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}
		for (iy=0; iy < nky; iy++) {
		    kiss_fft_stride(xinvs,(const kiss_fft_cpx *) (pp[0]+iy), 
				    (kiss_fft_cpx *) ctrace,nky);
		    for (ix=0; ix<nkx; ix++) {
			pp[ix][iy] = ctrace[ix];
		    }
		}

		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = cexpf(-csqrtf(w2*vt[ix][iy][iz])*dz);
			pp[ix][iy] = q[ix][iy][iz] + pp[ix][iy]*cshift;
		    }
		    for (iy=ny; iy<nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}
		for (ix=nx; ix<nkx; ix++) {
		    for (iy=0; iy<nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}
	    }

	    for (ix=0; ix<nx; ix++) {
		for (iy=0; iy < ny; iy++) {
		    cp[ix][iy][iw] = pp[ix][iy];
		}
	    }
	} else { /* migration */
    
	    for (ix=0; ix<nx; ix++) {
		for (iy=0; iy<ny; iy++) {
		    pp[ix][iy] = cp[ix][iy][iw];
		}
		for (iy=ny; iy < nky; iy++) {
		    pp[ix][iy] = 0.;
		}
	    }
	    for (ix=nx; ix<nkx; ix++) {
		for (iy=0; iy<nky; iy++) {
		    pp[ix][iy] = 0.;
		}
	    }

	    /* loop over migrated depths z */
	    for (iz=0; iz<nz-1; iz++) {
		/* accumulate image (summed over frequency) */
		for (ix=0; ix<nx; ix++) { 
		    for (iy=0; iy<ny; iy++) {
			q[ix][iy][iz] += crealf(pp[ix][iy]);
		    }
		}
		
		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = conjf(
			    cexpf(-0.5*csqrtf(w2*vt[ix][iy][iz])*dz));
			pp[ix][iy] *= cshift;
		    }
		    for (iy=ny; iy < nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}
		for (ix=nx; ix<nkx; ix++) {
		    for (iy=0; iy < nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}

		for (iy=0; iy < nky; iy++) {
		    kiss_fft_stride(xforw,(const kiss_fft_cpx *) (pp[0]+iy), 
				    (kiss_fft_cpx *) ctrace,nky);
		    for (ix=0; ix<nkx; ix++) {
			pp[ix][iy] = ctrace[ix];
		    }
		}
		for (ix=0; ix < nkx; ix++) {
		    kiss_fft(yforw,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}

		for (ix=0; ix<nkx; ix++) {
		    jx = (ix < nkx/2)? ix + nkx/2: ix - nkx/2;
		    kx = fkx + jx*dkx;
		    kx *= kx;

		    for (iy=0; iy<nky; iy++) {
			jy = (iy < nky/2)? iy + nky/2: iy - nky/2;
			ky = fky + jy*dky;
			ky = ky*ky + kx;

			cshift = conjf(cexpf((0.5*(csqrtf(w2*v[iz])+
						   csqrtf(w2*v[iz+1]))-
					      csqrtf(w2*v[iz]+ky))*dz));
			pp[ix][iy] *= cshift/(nkx*nky);
			/* Fourier scaling included */
		    }
		}

		for (ix=0; ix < nkx; ix++) {
		    kiss_fft(yinvs,(const kiss_fft_cpx *) pp[ix], 
			     (kiss_fft_cpx *) pp[ix]);
		}
		for (iy=0; iy < nky; iy++) {
		    kiss_fft_stride(xinvs,(const kiss_fft_cpx *) (pp[0]+iy), 
				    (kiss_fft_cpx *) ctrace,nky);
		    for (ix=0; ix<nkx; ix++) {
			pp[ix][iy] = ctrace[ix];
		    }
		}

		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = conjf(
			    cexpf(-0.5*csqrtf(w2*vt[ix][iy][iz+1])*dz));
			pp[ix][iy] *= cshift;
		    }
		    for (iy=ny; iy<nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}
		for (ix=nx; ix<nkx; ix++) {
		    for (iy=0; iy<nky; iy++) {
			pp[ix][iy] = 0.;
		    }
		}
	    }

	    for (ix=0; ix<nx; ix++) { 
		for (iy=0; iy<ny; iy++) {
		    q[ix][iy][nz-1] += crealf(pp[ix][iy]);
		}
	    }
	}
    } /* iw */
    
    free (*pp);
    free (pp);
    free (ctrace);
    free (xforw);
    free (xinvs);
    free (yforw);
    free (yinvs);
}

/* 	$Id: split1.c 790 2004-09-10 18:51:51Z fomels $	 */
