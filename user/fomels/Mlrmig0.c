/* Fast zero-offset time migration */
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

#include "fft2.h"

int main(int argc, char* argv[])
{
    int iz,im,ik,ix,i,j;     /* index variables */
    int nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2;
    float c;

    sf_complex *cwavem;
    float **wave, *curr;

    sf_complex **lt, **rt;
    sf_file left, right, inp, mig;

    sf_init(argc,argv);
    inp = sf_input("in");
    mig = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");

    nk = fft2_init(false,1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_complexalloc2(nzx,m2);
    rt = sf_complexalloc2(m2,nk);

    sf_complexread(lt[0],nzx*m2,left);
    sf_complexread(rt[0],m2*nk,right);

    wave = sf_floatalloc2(nzx2,m2);
    cwavem = sf_complexalloc(nk);
    curr = sf_floatalloc(nz);

    ifft2_allocate(cwavem);

    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	for (ik = 0; ik < nk; ik++) {
	    cwavem[ik] = rt[ik][im];
	}
	
	ifft2(wave[im],cwavem);
    }

    for (ix = 0; ix < nx; ix++) {
#ifdef _OPENMP
#pragma omp parallel for private(iz,i,j,c,im)
#endif
	for (iz=0; iz < nz; iz++) {
	    i = iz+ix*nz;  /* original grid */
	    j = iz+ix*nz2; /* padded grid */

	    c = 0.0f;
		
	    for (im = 0; im < m2; im++) {
		c += crealf(lt[im][i])*wave[im][j];
	    }

	    curr[iz] = c;
	}
	    	
	sf_floatwrite(curr,nz,mig);
    }
    
    exit (0);
}
