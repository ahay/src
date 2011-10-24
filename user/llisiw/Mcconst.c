/* Gaussian beam and exact complex eikonal for constant velocity medium */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
/*
NOTE: the current code only works in 2D and source on surface.
*/
    int  n1, n2, iz, ix;
    float  o1, o2, d1, d2, source, s, v0, angle;
    double sine, cosine, z, x, real, imag, amp, pha;
    double pi=3.14159265358979323846264338327950288;
    char *what;
    sf_complex *output;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* NOTE: the input only provides dimensional information. */
    if (!sf_histint(in,"n1",&n1)) sf_warning("No n1 in input.");
    if (!sf_histint(in,"n2",&n2)) n2=1;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_histfloat(in,"d1",&d1)) sf_warning("No d1 in input.");
    if (!sf_histfloat(in,"d2",&d2)) d2 = d1;

    if (!sf_getfloat("source",&source)) source=o2;
    /* real source point (on surface) */

    if (!sf_getfloat("s",&s)) s=0.;
    /* complex source shift */

    if (!sf_getfloat("v0",&v0)) v0=1.;
    /* constant velocity background */

    if (!sf_getfloat("angle",&angle)) angle=45.;
    /* rotation angle (counter-clock wise with respect to vertically downward) */ 
 
    if (NULL == (what = sf_getstring("what"))) what="exact";
    /* what to compute (default exact solution) */
   
    output = sf_complexalloc(n1*n2);
    sf_settype(out,SF_COMPLEX);

    /* convert angle to sin and cos */
    sine = sin(angle/180.*pi);
    cosine = sqrt(1.-sine*sine);

    switch (what[0]) {
	case 'e': /* exact solution */
	    for (ix=0; ix < n2; ix++) {
		for (iz=0; iz < n1; iz++) {
		    
		    /* coordinate rotation and shift */
		    z = (o1+iz*d1)*cosine+(o2+ix*d2-source)*sine;
		    x = -(o1+iz*d1)*sine+(o2+ix*d2-source)*cosine;
		    
		    /* mirror */
		    if (z < 0.)
			z = -z;

		    /* real and imaginary traveltime in reference coordinate */
		    real = z*z-s*s+x*x;
		    imag = -2.*z*s;
		    
		    amp = sqrt(sqrt(real*real+imag*imag));
		    pha = 0.;
		    
		    /* choose correct branch of square-root */
		    if (real > 0. && imag > 0.)
			pha = atan(imag/real)/2.;
		    if (real < 0. && imag > 0.)
			pha = (atan(imag/real)+pi)/2.;
		    if (real < 0. && imag < 0.)
			pha = (atan(imag/real)-pi)/2.;
		    if (real > 0. && imag < 0.)
			pha = atan(imag/real)/2.;
		    if (real >= 0. && imag == 0.)
			pha = 0.;
		    if (real == 0. && imag > 0.)
			pha = pi/4.;
		    if (real < 0. && imag == 0.)
			pha = -pi/2.;
		    if (real == 0. && imag < 0.)
			pha = -pi/4.;
		    
		    output[iz+ix*n1] = sf_cmplx((float)(amp*cos(pha)/v0),(float)((amp*sin(pha)+s)/v0));
		}
	    }
	    
	    sf_complexwrite(output,n1*n2,out);
	    break;

	case 'b': /* beam */
	    for (ix=0; ix < n2; ix++) {
		for (iz=0; iz < n1; iz++) {

		    /* coordinate rotation and shift */
		    z = (o1+iz*d1)*cosine+(o2+ix*d2-source)*sine;
		    x = -(o1+iz*d1)*sine+(o2+ix*d2-source)*cosine;

		    /* mirroe */
		    if (z < 0.)
			z = -z;

		    output[iz+ix*n1] = sf_cmplx((float)((z+z*x*x/(2.*(z*z+s*s)))/v0),(float)((s*x*x/(2.*(z*z+s*s)))/v0));
		}
	    }
	    
	    sf_complexwrite(output,n1*n2,out);
	    break;
    }

    exit(0);
}
