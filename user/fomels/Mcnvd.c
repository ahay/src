/* Residual from convection filtering. */
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
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int order, n1, n2, n12, i1, i2, iw, nw;
    float **data, ***conv, **resd, dif;
    sf_file in, cnv, out;

    sf_init(argc, argv);
    in = sf_input("in");
    cnv = sf_input("cnv");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");
    n12 = n1*n2;
    
    if (!sf_histint(cnv,"n2",&nw)) sf_error("Need n2= in cnv");
    order = nw/2;
  
    data = sf_floatalloc2(n1,n2);
    resd = sf_floatalloc2(n1,n2);
    conv = sf_floatalloc3(n1,nw,n2);

    sf_floatread(data[0],n12,in);
    sf_floatread(conv[0][0],n12*nw,cnv);

    for (i2=0; i2 < n2-1; i2++) {
	resd[i2][0] = 0.0f;
	for (i1=1; i1 < n1-1; i1++) {
	    dif = data[i2+1][i1] - data[i2][i1];
	    resd[i2][i1] = -dif;
	    for (iw=1; iw <= order; iw++) {
		resd[i2][i1] += conv[i2][2*iw-2][i1]*(dif - data[i2+1][i1+iw] + data[i2][i1-iw]);
		resd[i2][i1] += conv[i2][2*iw-1][i1]*(dif - data[i2+1][i1-iw] + data[i2][i1+iw]);
	    }
	}
	resd[i2][n1-1] = 0.0f;
    }
    for (i1=0; i1 < n1; i1++) {
	resd[n2-1][i1] = 0.0f;
    }

    sf_floatwrite(resd[0],n12,out);
  
    exit(0);
}
  

  
