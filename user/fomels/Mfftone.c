/* Test 1-D Fourier transform. */
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

#include "fft1.h"

int main(int argc, char* argv[])
{   
    bool inv;
    int nz,nz2,nk,ix,n2;
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
		
	if (!sf_getint("n1",&nz)) sf_error("Need n1=");	
	/* dimension (for inv=y) */
	sf_putint(space,"n1",nz);
    } else {
	space = sf_input("in");
	freq = sf_output("out");
		
	if (SF_FLOAT != sf_gettype(space)) sf_error("Need float input");
	sf_settype(freq,SF_COMPLEX);
		
	if (!sf_histint(space,"n1",&nz)) sf_error("No n1= in input");
    }
	
    nk = fft1_init(nz,&nz2);

    f = sf_floatalloc(nz2);
    c = sf_complexalloc(nk);
	
    if (inv) {
	if (!sf_histint(freq,"n1",&n2) || n2 != nk) 
	    sf_error("Need n1=%d in input",nk);
	ifft1_allocate(c);
    } else {
	sf_putint(freq,"n1",nk);
    }
		
    if (inv) {
	sf_complexread(c,nk,freq);
		
	ifft1(f,c);
		
	sf_floatwrite(f,nz,space);
    } else {
	for (ix=nz; ix < nz2; ix++) {
	    f[ix]=0.;
	}
		
	sf_floatread(f,nz,space);
		
	fft1(f,c);
		
	sf_complexwrite(c,nk,freq);
    }
	
    exit(0);
}
