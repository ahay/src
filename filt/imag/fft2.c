/* 2-D FFT encapsulated */
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

#include <rsf.h>
/*^*/

#include "fft2.h"

static kiss_fft_cfg xforw, xinvs, yforw, yinvs;
static int nx, ny;
static float complex *ctrace;

void fft2_init(int ny1, int nx1 /* in-line, cross-line */)
/*< initialize >*/
{
    nx = nx1;
    ny = ny1;

    xforw = kiss_fft_alloc(nx,0,NULL,NULL);
    xinvs = kiss_fft_alloc(nx,1,NULL,NULL);
    yforw = kiss_fft_alloc(ny,0,NULL,NULL);
    yinvs = kiss_fft_alloc(ny,1,NULL,NULL);

    ctrace = sf_complexalloc(nx);

    if (NULL == xforw || NULL == xinvs || NULL == yforw || NULL == yinvs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);
}

void fft2_close(void)
/*< Free allocated storage >*/
{
    free (ctrace);
    free (xforw);
    free (xinvs);
    free (yforw);
    free (yinvs);
}

void fft2(bool inv           /* inverse/forward flag */, 
	  complex float **pp /* [nx][ny] */) 
/*< Apply 2-D FFT >*/
{
    int ix, iy;
    
    if (inv) {
	for (ix=0; ix < nx; ix++) {
	    kiss_fft(yinvs,
		     (const kiss_fft_cpx *) pp[ix], 
		     (      kiss_fft_cpx *) pp[ix]);
	}
	for (iy=0; iy < ny; iy++) {
	    kiss_fft_stride(xinvs,
			    (const kiss_fft_cpx *) (pp[0]+iy), 
			    (      kiss_fft_cpx *) ctrace,ny);
	    for (ix=0; ix<nx; ix++) {
		pp[ix][iy] = ctrace[ix]/(nx*ny);
	    }
	}
    } else {
	for (iy=0; iy < ny; iy++) {
	    kiss_fft_stride(xforw,
			    (const kiss_fft_cpx *) (pp[0]+iy), 
			    (      kiss_fft_cpx *) ctrace,ny);
	    for (ix=0; ix<nx; ix++) {
		pp[ix][iy] = ctrace[ix];
	    }
	}
	for (ix=0; ix < nx; ix++) {
	    kiss_fft(yforw,
		     (const kiss_fft_cpx *) pp[ix], 
		     (      kiss_fft_cpx *) pp[ix]);
	}
    }
}
