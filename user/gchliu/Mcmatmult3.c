/* Multiplication of two complex matrices for 3D data.
I(n1,n2,f)*I(n2,n3,f)=D(n1,n3,f)
*/
/*
  Copyright (C) 2005 University of Texas at Austin

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
    int n1, n2, n3,n4,n5, i1, i2, i3,i4;
    float d2, o2;
    sf_complex **x, **y, **a;
    sf_file in=NULL, out=NULL, mat=NULL;

    sf_init(argc,argv);
    in  = sf_input ("in");
    out = sf_output("out");
    mat = sf_input ("mat");

    if (SF_COMPLEX != sf_gettype(in) ||
	SF_COMPLEX != sf_gettype(mat)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n4 = sf_leftsize(in,2);
    n5 = sf_leftsize(mat,2);
    if (n4 != n5) sf_error("Need same left dimention");
    x = sf_complexalloc2(n1,n2);

    if (!sf_histint(mat,"n1",&n3) || n3 != n2) sf_error("Need n1=%d in mat",n2);
    if (!sf_histint(mat,"n2",&n3)) sf_error("No n2= in mat");

    if (!sf_histfloat(mat,"d2",&d2)) d2=1.;
    if (!sf_histfloat(mat,"o2",&o2)) o2=0.;
    a = sf_complexalloc2(n2,n3);

    sf_putint(out,"n2",n3);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o2",o2);

    y = sf_complexalloc2(n1,n3);

    for (i4 = 0; i4 < n4; i4++) {
        sf_warning("slice %d of %d;",i4+1,n4);
        sf_complexread(x[0],n1*n2,in);
        sf_complexread(a[0],n2*n3,mat);
    	for (i1 = 0; i1 < n1; i1++) {
    	for (i3 = 0; i3 < n3; i3++) {
	         y[i3][i1] = sf_cmplx(0.,0.);
	    }
	    }
        for (i3 = 0; i3 < n3; i3++) {
	        for (i2 = 0; i2 < n2; i2++) {
	        for (i1 = 0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        y[i3][i1] += a[i3][i2] * x[i2][i1];
#else
		        y[i3][i1] = sf_cadd(y[i1],sf_cmul(a[i3][i2],x[i2][i1]));
#endif
	        }
	        }
	    }
	    sf_complexwrite(y[0],n1*n3,out);
    }
    exit(0);
}

