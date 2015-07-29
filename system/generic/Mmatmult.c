/* Simple matrix multiplication */
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
    bool adj;
    int n1, n2, i1, i2;
    float *x=NULL, *y=NULL, **a=NULL;
    sf_file in=NULL, out=NULL, mat=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    mat = sf_input("mat");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(mat)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

    if (!sf_getbool("adj",&adj)) adj=false;

    if (adj) {
	if (!sf_histint(mat,"n2",&n2) || n2 != n1) sf_error("Need n2=%d in mat",n1);
	if (!sf_histint(mat,"n1",&n2)) sf_error("No n1= in mat");
	a = sf_floatalloc2(n2,n1);
    } else {
	if (!sf_histint(mat,"n1",&n2) || n2 != n1) sf_error("Need n1=%d in mat",n1);
	if (!sf_histint(mat,"n2",&n2)) sf_error("No n2= in mat");
	a = sf_floatalloc2(n1,n2);
    }

    sf_putint(out,"n1",n2);

    x = sf_floatalloc(n1);
    y = sf_floatalloc(n2);

    sf_floatread(x,n1,in);
    sf_floatread(a[0],n1*n2,mat);

    for (i2 = 0; i2 < n2; i2++) {
	y[i2] = 0.;
	for (i1 = 0; i1 < n1; i1++) {
	    if (adj) {
		y[i2] += a[i1][i2] * x[i1];
	    } else {
		y[i2] += a[i2][i1] * x[i1];
	    }
	}
    }

    sf_floatwrite(y,n2,out);

    exit(0);
}
