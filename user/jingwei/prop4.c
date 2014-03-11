/* First nsps(+) then pspi(-), half time step */
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
#include "cfft2nsps.h"

/*******************************************************/

int prop4(sf_complex *input, sf_complex *output, sf_complex *lt, sf_complex *rt, int nz, int nx, int nkzx, int m2)
/*< First nsps(+) then pspi(-) >*/
{
    int iz, ix, im, ik, i, j;
    int nz2, nx2, nk, nzx, nzx2;
    int pad1 = 1;
    sf_complex **wave, **wave2, *curr, *currm, *cwave, *cwavem, c;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    if (nk!=nkzx) sf_error("nk discrepancy!");
    
    nzx = nz*nx;
    nzx2 = nz2*nx2;

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

    /* initialization */
    for (ix = 0; ix < nx2; ix++) {
	for (iz=0; iz < nz2; iz++) {
            i = iz+ix*nz;
	    j = iz+ix*nz2;
	    if (ix<nx && iz<nz)
		curr[j] = input[i];
	    else 
		curr[j] = sf_cmplx(0.,0.);
	}
    }


        /* nsps(+) */

        /* matrix multiplication */
	for (im = 0; im < m2; im++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		    currm[j] = lt[im*nzx+i]*curr[j];
#else
		    currm[j] = sf_cmul(lt[im*nzx+i], curr[j]);
#endif
		}
	    }
	    cfft2(currm,wave[im]);
	}
	
	for (ik = 0; ik < nk; ik++) {
	    c = sf_cmplx(0.,0.);
	    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		c += wave[im][ik]*rt[ik*m2+im];
#else
		c += sf_cmul(wave[im][ik],rt[ik*m2+im]);
#endif
	    }
	    cwave[ik] = c;
	}
	

        /* pspi(-) */
         
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	        cwavem[ik] = cwave[ik]*conjf(rt[ik*m2+im]);
#else
		cwavem[ik] = sf_cmul(cwave[ik],conjf(rt[ik*m2+im]));
#endif
	    }
	    icfft2(wave2[im],cwavem);
	}
	
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += conjf(lt[im*nzx+i])*wave2[im][j];
#else
		    c += sf_cmul(conjf(lt[im*nzx+i]), wave2[im][j]);
#endif
		}
		curr[j] = c;
	    }
	}       


    /* output final result*/
    for (ix = 0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
            i = iz+ix*nz;
	    j = iz+ix*nz2;
	    output[i] = curr[j];
	}
    }
    
    cfft2_finalize();
    return 0;
}
