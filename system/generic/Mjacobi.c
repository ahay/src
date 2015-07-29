/* Find eigenvalues of a symmetric matrix by Jacobi's iteration. */
/*
  Copyright (C) 2007 The University of Texas at Austin

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
#include "jacobi.h"

int main(int argc, char* argv[])
{
    int j, k, n, n2, i3, n3, iter, niter;
    float **a, *e, **v, s2;
    sf_file mat, val, eig;

    sf_init(argc,argv);
    mat = sf_input("in");
    val = sf_output("out");

    if (SF_FLOAT != sf_gettype(mat)) sf_error("Need float input");
    if (!sf_histint(mat,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histint(mat,"n2",&n2) || n2 != n) sf_error("Need n1=n2 in input");
    n3 = sf_leftsize(mat,2);

    sf_putint(val,"n2",1);

    if (!sf_getint("niter",&niter)) niter=10;

    a = sf_floatalloc2(n,n);
    e = sf_floatalloc(n);

    if (NULL != sf_getstring("eig")) {
	eig = sf_output("eig"); /* eigenvectors */
	v = sf_floatalloc2(n,n);
	for (j=0; j < n; j++) {
	    for (k=0; k < n; k++) {
		v[j][k] = (j==k)? 1.0:0.0;
	    }
	}
    } else {
	eig = NULL;
	v = NULL;
    }

    jacobi_init(n);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(a[0],n*n,mat);
	
	for (iter=0; iter < niter; iter++) {
	    s2 = 0.;
	    for (j=0; j < n-1; j++) {
		for (k=j+1; k < n; k++) {
		    s2 += jacobi(a,j,k,v);
		}
	    }
	    sf_warning("iter=%d s2=%g",iter+1,s2);
	}

	for (j=0; j < n; j++) {
	    e[j]=a[j][j];
	}

	sf_floatwrite(e,n, val);
	if (NULL != v) 	sf_floatwrite(v[0],n*n, eig);
    }

    exit(0);
}
