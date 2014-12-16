/* Output interpolation filter. 

See also: inttest1.
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
#include <string.h>
 
#include <rsf.h>

#include "interp_cube.h"
#include "interp_sinc.h"
#include "interp_mom.h"

int main(int argc, char* argv[])
{
    int nw;
    float *w, kai, x;
    char *intp;
    sf_interpolator interp=NULL;
    sf_file filt;

    sf_init (argc,argv);
    filt = sf_output("out");
    sf_setformat(filt,"native_float");

    intp = sf_getstring("interp");
    /* interpolation (lagrange,cubic,kaiser,lanczos,cosine,welch,spline,mom) */
    if (NULL == intp) sf_error("Need interp=");

    if (!sf_getint("nw",&nw)) sf_error("Need nw=");
    /* interpolator size */

    sf_putint(filt,"n1",nw);
    w = sf_floatalloc(nw);

    if (!sf_getfloat("x",&x)) sf_error("Need x=");
    /* interpolation shift */

    switch(intp[0]) {
	case 'l':
	    if (!strncmp("lag",intp,3)) { /* Lagrange */
		interp = sf_lg_int;
	    } else if (!strncmp("lan",intp,3)) { /* Lanczos */
		sinc_init('l', 0.);
		interp = sinc_int;
	    } else {
		sf_error("%s interpolator is not implemented",intp);
	    }
	    break;
	case 'c':
	    if (!strncmp("cub",intp,3)) { /* Cubic convolution */
		interp = sf_lg_int;
	    } else if (!strncmp("cos",intp,3)) { /* Cosine */
		sinc_init('c', 0.);
		interp = sinc_int;
	    } else {
		sf_error("%s interpolator is not implemented",intp);
	    }
	    break;
	case 'k':
	    if (!sf_getfloat("kai",&kai)) kai=4.0;
	    /* Kaiser window parameter */
	    sinc_init('k', kai);
	    interp = sinc_int;
	    break;
	case 'w':
	    sinc_init('w', 0.);
	    interp = sinc_int;
	    break;
	case 's':
	    interp = sf_spline_int;
	    break;
	case 'm':
	    interp = mom_int;
	    break;
	case 'h':
	    x-=0.21; // optimal shift
	    interp = sf_lin_int;
	    break;
	default:
	    sf_error("%s interpolator is not implemented",intp);
	    break;
    }

    interp(x,nw,w);
    sf_floatwrite(w,nw,filt);

    exit(0);
}



