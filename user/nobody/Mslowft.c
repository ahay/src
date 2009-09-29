/* Slow FT transform on the first axis.

Input and output are complex data. The input is padded by factor pad.
*/
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

int main (int argc, char **argv)
{
    int nx, n1, n2;	     /* dimensions */
    int ix, i2, ik;       /* loop counters 	*/
    int nk;		     /* number of wavenumbers */	
    int npad;                /* padding */

    float dx;		     /* space sampling interval */
    float dk;	             /* wavenumber sampling interval */
    float x0;                /* staring space */
    float k0;                /* starting wavenumber */
    float wt;                /* Fourier scaling */
    float x, k;

    sf_complex *ck, *cx;     /* frequency-wavenumber */
    sf_complex shift;        /* phase shift */

    bool inv;                /* forward or inverse */
    bool sym;                /* symmetric scaling */
    int sign;                /* transform sign */

    sf_file in, out;

    sf_init(argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error ("Need complex input");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* if y, perform inverse transform */

    if (!sf_getbool("sym",&sym)) sym=false;
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */

    if (!sf_getint("sign",&sign)) sign = inv? 1: 0;
    /* transform sign (0 or 1) */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    
    n2 = sf_leftsize(in,1);
    
    if (inv) { 
	nk = n1;
	if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");
	if (!sf_histint  (in,"m1",&nx)) nx=nk;
	if (!sf_histfloat(in,"c1",&x0)) x0 = 0.; 
	dx = 1./(nk*dk);

	sf_putint (out,"n1",nx);
	sf_putfloat (out,"d1",dx);
	sf_putfloat (out,"o1",x0);
    } else { 
	nx = n1;
	if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&x0)) x0 = 0.;

	sf_putint(out,"m1",nx);
	sf_putfloat(out,"c1",x0);

	if (!sf_getint("pad",&npad)) npad=2;
	/* padding factor */

	/* determine wavenumber sampling */
	nk = nx*npad;
	dk = 1./(nk*dx);
	k0 = -0.5/dx;

	sf_putint (out,"n1",nk);
	sf_putfloat (out,"d1",dk);
	sf_putfloat (out,"o1",k0);
    }
    

    ck = sf_complexalloc(nk);
    cx = sf_complexalloc(nx);

    /* FFT scaling */
    wt = sym? 1./sqrtf((float) nk): 1./nk;
    
    for (i2=0; i2<n2; i2++) {
	if (inv) {
	    sf_complexread(ck,nk,in);	    
	} else {
	    sf_complexread(cx,nx,in);
	}

	for (ik=0; ik<nk; ik++) {
	    k = k0+ik*dk;
	    for (ix=0; ix<nx; ix++) {
		x = x0+ix*dx;
		x = 2.*SF_PI*x*k;
		shift = sf_cmplx(cosf(x),sinf(x));

#ifdef SF_HAS_COMPLEX_H
		if (inv) {
		    cx[ix] = ck[ik]*shift*wt;
		} else {
		    ck[ix] = cx[ik]*shift*wt;
		}
#else
		if (inv) {
		    cx[ix] = sf_cmul(ck[ik],sf_crmul(shift,wt));
		} else {
		    ck[ix] = sf_cmul(cx[ik],sf_crmul(shift,wt));
		}
#endif
	    }
	}

	if (inv) {      
	    sf_complexwrite(cx,nx,out);
	} else {
	    sf_complexwrite(ck,nk,out);
	}
    }

    exit (0);
}

/* 	$Id: Mfft3.c 1151 2005-05-25 11:27:14Z fomels $	 */
