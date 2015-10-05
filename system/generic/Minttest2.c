/* Interpolation from a regular grid in 2-D. */
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
#include "shprefilter.h"

int main(int argc, char* argv[])
{
    int m[2], n, n3, nd1, nd2, nd, nw, i3, two;
    float **mm, **coord, *z, *tmp=NULL;
    float o1,o2, oo1,oo2, d1,d2, dd1,dd2, kai;
    char *intp;
    sf_interpolator interp;
    sf_bands spl1=NULL, spl2=NULL;
    sf_file in, out, crd;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    crd = sf_input("coord");

    if (!sf_histint(in,"n1",m))   sf_error("No n1= in input");
    if (!sf_histint(in,"n2",m+1)) sf_error("No n2= in input");
    n = m[0]*m[1];
    n3 = sf_leftsize(in,2);

    if (!sf_histint(crd,"n1",&two) || 2 != two) 
	sf_error("Need n1=2 in coord");
    if (!sf_histint(crd,"n2",&nd1)) sf_error("No n2= in coord");
    if (!sf_histint(crd,"n3",&nd2)) nd2=1;
    nd = nd1*nd2;
    sf_putint(out,"n1",nd1);
    sf_putint(out,"n2",nd2);

    if (!sf_histfloat(in,"d1",&d1))   sf_error("No d1= in input");
    if (!sf_histfloat(crd,"d2",&dd1)) dd1=d1;
    sf_putfloat(out,"d1",dd1);

    if (!sf_histfloat(in,"d2",&d2))   sf_error("No d2= in input");
    if (!sf_histfloat(crd,"d3",&dd2)) dd2=d2;
    sf_putfloat(out,"d2",dd2);

    if (!sf_histfloat(in,"o1",&o1))   sf_error("No o1= in input");
    if (!sf_histfloat(crd,"o2",&oo1)) oo1=o1;
    sf_putfloat(out,"o1",oo1);

    if (!sf_histfloat(in,"o2",&o2))   sf_error("No o2= in input");
    if (!sf_histfloat(crd,"o3",&oo2)) oo2=o2;
    sf_putfloat(out,"o2",oo2);

    intp = sf_getstring("interp");
    /* interpolation (lagrange,cubic,kaiser,lanczos,cosine,welch,spline) */
    if (NULL == intp) sf_error("Need interp=");

    if (!sf_getint("nw",&nw)) sf_error("Need nw=");
    /* interpolator size */

    coord = sf_floatalloc2(2,nd);
    sf_floatread(coord[0],nd*2,crd);

    switch(intp[0]) {
	case 'l':
	    if (!strncmp("lag",intp,3)) { /* Lagrange */
		interp = sf_lg_int;
	    } else if (!strncmp("lan",intp,3)) { /* Lanczos */
		sinc_init('l', 0.);
		interp = sinc_int;
	    } else {
		interp = NULL;
		sf_error("%s interpolator is not implemented",intp);
	    }
	    break;
	case 'c':
	    if (!strncmp("cub",intp,3)) { /* Cubic convolution */
		interp = cube_int;
	    } else if (!strncmp("cos",intp,3)) { /* Cosine */
		sinc_init('c', 0.);
		interp = sinc_int;
	    } else {
		interp = NULL;
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
	    spl1 = sf_spline_init(nw,m[0]);
	    spl2 = sf_spline_init(nw,m[1]);
	    tmp = sf_floatalloc(m[1]);
	    break;
	case 'h':
	    interp = sf_lin_int;
	    break;
	default:
	    interp = NULL;
	    sf_error("%s interpolator is not implemented",intp);
	    break;
    }

    if (intp[0] == 'h') sf_int2sh_init (coord, o1, o2, d1, d2, m[0], m[1], interp, nw, nd);
    else sf_int2_init (coord, o1, o2, d1, d2, m[0], m[1], interp, nw, nd);

    z = sf_floatalloc(nd);
    mm = sf_floatalloc2(m[0],m[1]);

    for (i3=0; i3 < n3; i3++) {
        sf_floatread (mm[0],n,in);
        if ('s' == intp[0]) sf_spline2(spl1,spl2,m[0],m[1],mm,tmp);
        if ('h' == intp[0]) shprefilter2d(m[0],m[1],mm); 

	sf_int2_lop (false,false,n,nd,mm[0],z);

	sf_floatwrite (z,nd,out);
    }

    exit(0);
}
