/* Test 3-D Fourier transform. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "fft3.h"

int main(int argc, char* argv[])
{   
    bool inv, cmplx;
    int nz,nx,ny,nz2,nx2,ny2,nk,ix,iy,n2,pad1;
    float *f;
    sf_complex *c; 
    sf_file space, freq;

    sf_init(argc,argv);
    
    if (!sf_getbool("inv",&inv)) inv=false;
    /* inverse flag */

    if (inv) {
	freq = sf_input("in");
	space = sf_output("out");

	if (SF_COMPLEX != sf_gettype(freq)) sf_error("Need complex input");
	sf_settype(space,SF_FLOAT);

	if (!sf_getint("n1",&nz)) sf_error("No n1= in input");
	if (!sf_getint("n2",&ny)) sf_error("No n2= in input");
	if (!sf_getint("n3",&nx)) sf_error("No n3= in input");

	sf_putint(space,"n1",nz);
	sf_putint(space,"n2",ny);
	sf_putint(space,"n3",nx);
    } else {
	space = sf_input("in");
	freq = sf_output("out");
    
	if (SF_FLOAT != sf_gettype(space)) sf_error("Need float input");
	sf_settype(freq,SF_COMPLEX);
	
	if (!sf_histint(space,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histint(space,"n2",&ny)) sf_error("No n2= in input");
	if (!sf_histint(space,"n3",&nx)) sf_error("No n3= in input");
    }

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    nk = fft3_init(cmplx,pad1,nz,ny,nx,&nz2,&ny2,&nx2);

    f = sf_floatalloc(nz2*ny2*nx2);
    c = sf_complexalloc(nk);

    if (inv) {
	if (!sf_histint(freq,"n1",&n2) || n2 != nk) sf_error("Need n1=%d in input",nk);
	ifft3_allocate(c);
    } else {
	sf_putint(freq,"n1",nk);
	sf_putint(freq,"n2",1);
	sf_putint(freq,"n3",1);
    }

    if (inv) {
	sf_complexread(c,nk,freq);

	ifft3(f,c);

	for (ix=0; ix < nx; ix++) {
	    for (iy=0; iy < ny; iy++) {
		sf_floatwrite(f+(ix*ny2+iy)*nz2,nz,space);
	    }
	}
    } else {
	for (ix=0; ix < nz2*ny2*nx2; ix++) {
	    f[ix]=0.;
	}
	
	for (ix=0; ix < nx; ix++) {
	    for (iy=0; iy < ny; iy++) {
		sf_floatread(f+(ix*ny2+iy)*nz2,nz,space);
	    }
	}

	fft3(f,c);
	
	sf_complexwrite(c,nk,freq);
    }

    exit(0);
}
