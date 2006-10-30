/* Kroneker product with square matrices */
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

#include "kron.c"

int main(int argc, char* argv[])
{
    bool adj, inv;
    int n, n2, niter;
    float *x, *y, *p, **a, **b, eps;
    sf_file in, out, mat1, mat2;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    mat1 = sf_input("mat1");
    mat2 = sf_input("mat2");
    
    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(mat1) ||
	SF_FLOAT != sf_gettype(mat2)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2) || n2 != n) sf_error("Need n2=%d",n);
    n2 = n*n;

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */
    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;

    a = sf_floatalloc2(n,n);
    b = sf_floatalloc2(n,n);

    x = sf_floatalloc(n2);
    y = sf_floatalloc(n2);

    kron_init(n,a,b);

    sf_floatread(x,n2,in);
    sf_floatread(a[0],n2,mat1);
    sf_floatread(b[0],n2,mat2);

    if (adj) {
	if (inv) {
	    p = sf_floatalloc(n2);

	    sf_conjgrad_init(n2,n2,n2,n2,eps,FLT_EPSILON,true,false);
	    sf_conjgrad(NULL,kron_lop,sf_copy_lop,p,y,x,niter);
	} else {
	    kron_lop(true,false,n2,n2,y,x);
	} 
    } else {
	kron_lop(false,false,n2,n2,x,y);
    }

    sf_floatwrite(y,n2,out);

    exit(0);
}
