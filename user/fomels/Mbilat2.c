/* 2-D bilateral filtering */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
    int n1, n2, n12, i3, n3, i1, i2, r1, r2, k1, k2, j1, j2, repeat, ir;
    float **dinp, **dout, **norm, a1, a2, a3, d1, d2, d0, d, delta1, delta2, delta3, gauss;
    sf_file inp, out;
    
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(inp,2);

    if (!sf_histfloat(inp,"d1",&d1)) d1=1.0f;
    if (!sf_histfloat(inp,"d2",&d2)) d2=1.0f;

    dinp = sf_floatalloc2(n1,n2);
    dout = sf_floatalloc2(n1,n2);
    norm = sf_floatalloc2(n1,n2);

    if (!sf_getint("r1",&r1)) r1=1; /* vertical smoothing radius */
    if (!sf_getint("r2",&r2)) r2=1; /* horizontal smoothing radius */

    if (!sf_getfloat("a1",&a1)) a1=0.0f; /* vertical attenuation factor */
    if (!sf_getfloat("a2",&a2)) a2=a1;   /* horizontal attenuation factor */
    if (!sf_getfloat("a3",&a3)) a3=0.0f; /* data attenuation factor */ 

    if (!sf_getint("repeat",&repeat)) repeat=1; /* repeat the operation */

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(dinp[0],n12,inp);

	for (ir=0; ir < repeat; ir++) {

	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    dout[i2][i1] = 0.0f;
		    norm[i2][i1] = 0.0f;
		    d0 = dinp[i2][i1];

		    for (k2=-r2; k2 <= r2; k2++) {
			j2 = i2+k2;
			if (j2 < 0 || j2 >= n2) continue;
			delta2 = (j2-i2)*d2;
			delta2 *= delta2;

			for (k1=-r1; k1 <= r1; k1++) {
			    j1 = i1+k1;
			    if (j1 < 0 || j1 >= n1) continue;
			    delta1 = (j1-i1)*d1;
			    delta1 *= delta1;

			    d = dinp[j2][j1];
			    delta3 = d-d0;
			    delta3 *= delta3;

			    gauss = expf(-a1*delta1-a2*delta2-a3*delta3);

			    dout[i2][i1] += d*gauss;
			    norm[i2][i1] += gauss;
			}
		    }

		    if (norm[i2][i1] > 0.0f) {
			dinp[i2][i1] = dout[i2][i1]/norm[i2][i1]; 
		    } else {
			dinp[i2][i1] = 0.0f;
		    }
		}
	    }

	}
	
	sf_floatwrite(dinp[0],n12,out);
    }

    exit(0);
}
