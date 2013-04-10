/* Tridiagonal matrix solver using chasing method
Ax=d-> LUx=d -> Ly=d -> Ux=y -> x
*/
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

#include "tridsolver.h"
/*Ax=d*/
/* A is a toeplitz matrix but not limitted to a symmetric tridiagonal matrix. */

int main (int argc, char* argv[])
{

    int i,n1,n2,n22,n;
    float *A, *x, *d;
    bool verb;
    sf_file in, rhs, out;

    sf_init (argc,argv);
    in = sf_input("in");
    rhs= sf_input("rhs");
    out = sf_output("out");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input matrix!");
    if(!sf_histint(in,"n2",&n2)) sf_error("No n2 in input matrix!");
    if(!sf_histint(rhs,"n2",&n22))sf_error("No n2 in input vector!");
    if(!sf_getbool("verb",&verb)) verb=false;
    if((n1!=n2) || (n1!=n22)) sf_error("Dimension mistake!");
    n=n1;

    x = sf_floatalloc(n);
    A=sf_floatalloc(n*n);
    d=sf_floatalloc(n);

    sf_floatread(A,n*n,in);
    sf_floatread(d,n,rhs);

    trid_init(n, A);
    trid_solve(d, x);
    
    if(verb)
    {
	for(i=0;i<n;i++)
	sf_warning("x[%d]=%f",i+1,x[i]);
    }
    
   sf_putint(out,"n2",n);
   sf_putint(out,"n1",1);
   sf_floatwrite(x,n,out);

   trid_close();
   free(A);
   free(d);
   free(x);
   exit(0);
}

