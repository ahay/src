/* Test Beam for constant velocity gradient */
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
    int  n1, n2,*m;
    float  o1, o2, d1, d2, source, w, p, v0, b, z0, x0, r0, z, x, r, ztemp, xtemp, xi, lg;
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

    if (!sf_getfloat("w",&w)) w=0.;
    /* beam width */

    if (!sf_getfloat("v0",&v0)) v0=1.;
    /* surface velocity */
    
    if (!sf_getfloat("b",&b)) b=0.;
    /* velocity gradient */
 
    if (!sf_getfloat("p",&p)) p=1/v0;;
    /* ray parameter */

    output = sf_complexalloc(n1*n2);
    m = sf_intalloc(n1*n2);

    sf_settype(out,SF_COMPLEX);
    sf_settype(mask,SF_INT); 

    x0 = source-sqrtf(1/p/p-v0*v0)/b;
    z0 = -v0/b;
    r0 = 1/b/p;

    sf_warning("x0=%g z0=%g r0=%g",x0,z0,r0);

    for (ix=0; ix < n2; ix++) {
	x = o2+ix*d2;
	for (iz=0; iz < n1; iz++) {
	    z = o1+iz*d1;

	    r = hypotf(x-x0,z-z0);
	    xtemp = (x-x0)*r0/r+x0;
	    ztemp = (z-z0)*r0/r+z0;

	    xi = b*b*((xtemp-source)*(xtemp-source)+(ztemp*ztemp))/(2*(v0+b*ztemp)*v0);
	    if (x == source && z==0)
	    {
		output[iz+ix*n1] = sf_cmplx(0.,0.);
		m[iz+ix*n1] = 0;
	    }
	    else {
		lg = logf(1+xi+sqrtf((2+xi)*xi));
		output[iz+ix*n1] = sf_cmplx(1/b*lg+w*(r-r0)*(r-r0),
					    w*(r-r0)*(r-r0)*b/(lg+SF_SIG(lg)*0.001));
		if (w*(r-r0)*(r-r0)*b/logf(1+xi+sqrtf((2+xi)*xi))+0.001 <= 0.0011)
		    m[iz+ix*n1] = 0;
		else
		    m[iz+ix*n1] = 1;
	    }
	}
    }
    sf_complexwrite(output,n1*n2,out);
    sf_intwrite(m,n1*n2,mask);
    
    exit(0);
}
