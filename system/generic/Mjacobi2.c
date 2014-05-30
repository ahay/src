/* Find eigenvalues of a general complex matrix by Jacobi-like iteration. */
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
    bool verb;
    int j, k, n, n2, i3, n3, iter, niter;
    sf_complex **a=NULL, *e=NULL;
    float s2;
    sf_file mat=NULL, val=NULL;

    sf_init(argc,argv);
    mat = sf_input("in");
    val = sf_output("out");

    if (SF_COMPLEX != sf_gettype(mat)) sf_error("Need complex input");
    if (!sf_histint(mat,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histint(mat,"n2",&n2) || n2 != n) sf_error("Need n1=n2 in input");
    n3 = sf_leftsize(mat,2);

    sf_putint(val,"n2",1);

    if (!sf_getint("niter",&niter)) niter=10;
    if (!sf_getbool("verb",&verb)) verb=false;

    a = sf_complexalloc2(n,n);
    e = sf_complexalloc(n);
    jacobi2_init(n,verb);

    for (i3=0; i3 < n3; i3++) {
	sf_complexread(a[0],n*n,mat);
	
	for (iter=0; iter < niter; iter++) {
	    s2 = 0.;
	    for (j=0; j < n; j++) {
		for (k=0; k < n; k++) {
		    s2 += jacobi2(a,n,j,k);
		}
	    }
	    sf_warning("iter=%d s2=%g",iter+1,s2);
	}

	for (j=0; j < n; j++) {
	    e[j]=a[j][j];
	}
	
	sf_complexwrite(e,n, val);
    }

    exit(0);
}
