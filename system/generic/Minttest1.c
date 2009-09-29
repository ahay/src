/* Testing forward interpolation in 1-D. */
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
    int n, n2, nd, nw, i2;
    float *mm=NULL, *coord=NULL, *z=NULL, o, oo, d, dd, kai;
    char *intp=NULL;
    sf_interpolator interp=NULL;
    sf_file in=NULL, out=NULL, crd=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    crd = sf_input("coord");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(crd,"n1",&nd)) sf_error("No n1= in coord");
    sf_putint(out,"n1",nd);

    if (!sf_histfloat(in,"d1",&d))   sf_error("No d1= in input");
    if (!sf_histfloat(crd,"d1",&dd)) sf_error("No d1= in coord");
    sf_putfloat(out,"d1",dd);

    if (!sf_histfloat(in,"o1",&o))   sf_error("No o1= in input");
    if (!sf_histfloat(crd,"o1",&oo)) sf_error("No o1= in coord");
    sf_putfloat(out,"o1",oo);

    intp = sf_getstring("interp");
    /* interpolation (lagrange,cubic,kaiser,lanczos,cosine,welch,spline,mom) */
    if (NULL == intp) sf_error("Need interp=");

    if (!sf_getint("nw",&nw)) sf_error("Need nw=");
    /* interpolator size */

    coord = sf_floatalloc(nd);
    sf_floatread(coord,nd,crd);

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
	    sf_prefilter_init (nw, n, 3*n);
	    interp = sf_spline_int;
	    break;
	case 'm':
	    sf_prefilter_init (-nw, n, 3*n);
	    interp = mom_int;
	    break;
	default:
	    sf_error("%s interpolator is not implemented",intp);
	    break;
    }

    sf_int1_init (coord, o, d, n, interp, nw, nd);

    z = sf_floatalloc(nd);
    mm = sf_floatalloc(n);
 
    for (i2=0; i2 < n2; i2++) {
        sf_floatread (mm,n,in);

        if ('s' == intp[0] || 'm' == intp[0]) 
	    sf_prefilter_apply (n,mm);
	
	sf_int1_lop (false,false,n,nd,mm,z);

	sf_floatwrite (z,nd,out);
    }

    exit(0);
}
