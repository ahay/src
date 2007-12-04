/* Find roots of a complex polynomial. */
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

#include "jacobi2.h"

int main(int argc, char* argv[])
{
    int j, k, n, m, i2, n2, iter, niter;
    sf_complex **a, *e;
    float s2;
    sf_file poly, root;

    sf_init(argc,argv);
    poly = sf_input("in");
    root = sf_output("out");

    if (SF_COMPLEX != sf_gettype(poly)) sf_error("Need complex input");
    if (!sf_histint(poly,"n1",&n)) sf_error("No n1= in input");
    n2 = sf_leftsize(poly,1);

    if (!sf_getint("niter",&niter)) niter=10;

    sf_putint(root,"n1",n-1);

    a = sf_complexalloc2(n,n);
    e = sf_complexalloc(n);
    jacobi2_init(n);

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(e,n,poly);

	for (m = n; m > 0; m--) {
	    if (cabsf(e[m-1]) > FLT_EPSILON) break;
	}
	m--;

	for (j=0; j < m; j++) {
	    for (k=0; k < m; k++) {
		a[j][k] = sf_cmplx(0.,0.);
	    }
	}
	for (j=0; j < m-1; j++) {
	    a[j][j+1]=sf_cmplx(1.,0.);
	}
	for (j=0; j < m; j++) {
#ifdef SF_HAS_COMPLEX_H
	    a[m-1][j]=-e[j]/e[m];
#else
	    a[m-1][j]=sf_cneg(sf_cdiv(e[j],e[m]));
#endif
	}

	for (iter=0; iter < niter; iter++) {
	    s2 = 0.;
	    for (j=0; j < m; j++) {
		for (k=0; k < m; k++) {
		    s2 += jacobi2(a,m,j,k);
		}
	    }
	    sf_warning("iter=%d s2=%g",iter+1,s2);
	}

	for (j=0; j < m; j++) {
	    e[j]=a[j][j];
	}
	for (j=m; j < n-1; j++) {
	    e[j]=sf_cmplx(0.,0.);
	}
	
	sf_complexwrite(e,n-1, root);
    }
    
    exit(0);
}
    

