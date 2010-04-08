/* Test Exact and Beam for constant velocity background */
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

int main(int argc, char* argv[])
{
    int  n1, n2, *m;
    float  o1, o2, d1, d2, source, s, v0, sin, cos, z, x, real, imag;
    char *what;
    int iz, ix;
    sf_complex *output;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("mask");

    if (!sf_histint(in,"n1",&n1)) sf_warning("No n1 in input.");
    if (!sf_histint(in,"n2",&n2)) n2=1;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_histfloat(in,"d1",&d1)) sf_warning("No d1 in input.");
    if (!sf_histfloat(in,"d2",&d2)) d2 = d1;

    if (!sf_getfloat("source",&source)) source=o2;
    /* real source point */

    if (!sf_getfloat("s",&s)) s=0.;
    /* imaginary source point */

    if (!sf_getfloat("v0",&v0)) v0=1.;
    /* constant velocity background */

    if (!sf_getfloat("sin",&sin)) sin=sqrtf(2)/2;
    /* angle */ 
 
    if (NULL==(what=sf_getstring("what"))) what="exact";
    /* what to compute */
   
    output = sf_complexalloc(n1*n2);
    m = sf_intalloc(n1*n2);

    sf_settype(out,SF_COMPLEX);
    sf_settype(mask,SF_INT);

    cos = sqrtf(1-sin*sin);

    switch (what[0]) {
	case 'e': /* exact solution */
	    for (ix=0; ix < n2; ix++) {
		for (iz=0; iz < n1; iz++) {
/*
		    x = 1/sqrtf(2)*(o2+ix*d2-source+o1+iz*d1);
		    z = 1/sqrtf(2)*(o1+iz*d1-o2-ix*d2+source);
*/
		    x = (o2+ix*d2-source)*cos+(o1+iz*d1)*sin;
		    z = (o1+iz*d1)*cos-(o2+ix*d2-source)*sin;
		    
		    real = x*x+z*z-s*s;
		    imag = -2*z*s;
		    
		    if (z > 0.)
		    {
			if ((imag/sqrtf(2*(hypotf(real,imag)+real))+s)/v0 >= 0.)
			    output[iz+ix*n1] = sf_cmplx(sqrtf((hypotf(real,imag)+real)/2)/v0,(imag/sqrtf(2*(hypotf(real,imag)+real))+s)/v0+0.001);
			else
			    output[iz+ix*n1] = output[iz+ix*n1-1];
		    }
		    if (z == 0.)
		    {
			if (real >= 0.)
			    output[iz+ix*n1] = sf_cmplx(sqrtf(real)/v0,s/v0);
			else
			    output[iz+ix*n1] = sf_cmplx(0.,(-sqrtf(-real)+s)/v0);
		    }
		    if (z < 0.)
		    {
			output[iz+ix*n1] = sf_cmplx(sqrtf((hypotf(real,imag)+real)/2)/v0,(-imag/sqrtf(2*(hypotf(real,imag)+real))+s)/v0+0.001);
		    }
		}
	    }
	    sf_complexwrite(output,n1*n2,out);
	    break;
	case 'b': /* beam */
	    for (ix=0; ix < n2; ix++) {
		for (iz=0; iz < n1; iz++) {
/*
		    x = 1/sqrtf(2)*(o2+ix*d2-source+o1+iz*d1);
		    z = 1/sqrtf(2)*(o1+iz*d1-o2-ix*d2+source);
*/
		    x = (o2+ix*d2-source)*cos+(o1+iz*d1)*sin;
		    z = (o1+iz*d1)*cos-(o2+ix*d2-source)*sin;

		    if (z == 0.)
		    {
			output[iz+ix*n1] = sf_cmplx(0.,x*x/2/s/v0);
			m[iz+ix*n1] = 0;
		    }
		    else {
			output[iz+ix*n1] = sf_cmplx((z+x*x*z/(2*(z*z+s*s)))/v0,(x*x*s/(2*(z*z+s*s)))/v0+0.001);
			if ((x*x*s/(2*(z*z+s*s)))/v0 <= 0.00001)
			    m[iz+ix*n1] = 0;
			else
			    m[iz+ix*n1] = 1;
		    }
		}
	    }
	    sf_complexwrite(output,n1*n2,out);
	    sf_intwrite(m,n1*n2,mask);
	    break;
    }

    exit(0);
}
