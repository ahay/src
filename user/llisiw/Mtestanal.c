/* Test Analytical for constant velocity background */
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
    int  n1, n2;
    float  o1, o2, d1, d2, source1, source2, s1, s2, s3, v0, real, imag;
    int ix, iy;
    sf_complex *output;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_warning("No n1 in input.");
    if (!sf_histint(in,"n2",&n2)) n2=1;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_histfloat(in,"d1",&d1)) sf_warning("No d1 in input.");
    if (!sf_histfloat(in,"d2",&d2)) d2 = d1;

    if (!sf_getfloat("source1",&source1)) source1=o1;
    if (!sf_getfloat("source2",&source2)) source2=o2;
    /* real source point */

    if (!sf_getfloat("s1",&s1)) s1=0.;
    if (!sf_getfloat("s2",&s2)) s2=0.;
    if (!sf_getfloat("s3",&s3)) s3=0.;
    /* imaginary source point */

    if (!sf_getfloat("v0",&v0)) v0=1.;
    /* constant velocity background */
    
    output = sf_complexalloc(n1*n2);

    sf_settype(out,SF_COMPLEX);

    if (s2==0. && s3==0.)
    {
	for (iy=0; iy < n2; iy++) {
	    for (ix=0; ix < n1; ix++) {
		if (hypotf(o1+ix*d1-source1,o2+iy*d2-source2) >= s1)
		    output[ix+iy*n1] = sf_cmplx(sqrtf((o1+ix*d1-source1)*(o1+ix*d1-source1)+(o2+iy*d2-source2)*(o2+iy*d2-source2)-s1*s1)/v0,s1/v0);
		else
		    output[ix+iy*n1] = sf_cmplx(0.,(-sqrtf(s1*s1-(o1+ix*d1-source1)*(o1+ix*d1-source1)-(o2+iy*d2-source2)*(o2+iy*d2-source2))+s1)/v0);
	    }
	}
    } else {
	for (iy=0; iy < n2; iy++) {
	    for (ix=0; ix < n1; ix++) {
		real = (o1+ix*d1-source1)*(o1+ix*d1-source1)+(o2+iy*d2-source2)*(o2+iy*d2-source2)-(s1*s1+s2*s2+s3*s3);
		imag = -2*(s2*(o1+ix*d1-source1)+s3*(o2+iy*d2-source2));
		if (imag != 0. || (imag == 0. && real >= 0.))
		{
		    output[ix+iy*n1] = 
			sf_cmplx(sqrtf((hypotf(real,imag)+real)/2)/v0,
				 (imag/sqrtf(2*(hypotf(real,imag)+real))+
				  sqrtf(s1*s1+s2*s2+s3*s3))/v0);
		} else {
		    output[ix+iy*n1] = sf_cmplx(0.,
						(-sqrtf(-real)+
						 sqrtf(s1*s1+s2*s2+s3*s3))/v0);
		}
	    }
	}
    }
    sf_complexwrite(output,n1*n2,out);
    exit(0);
}
