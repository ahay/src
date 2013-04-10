/* Symmetric positive definite matrix equation solver using square root method (cholesky decomposition method)
Ax=d-> LL'x=d -> Ly=d -> L'x=y -> x
Input is the down triangle of A.
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
#include "cholesky.h"
/*Ax=d*/
/* A is a symmetrix positive matrix */

int main (int argc, char* argv[])
{

    int i,n1,n11,n;
    float *A, *d;
    bool verb;
    sf_file in, rhs, out;

    sf_init (argc,argv);
    in = sf_input("in");
    rhs= sf_input("rhs");
    out = sf_output("out");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

   if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input matrix!");
   if(!sf_histint(rhs,"n1",&n11)) sf_error("No n1 in input vector!");
   if(!sf_getbool("verb",&verb)) verb=false;
   if(!sf_getint("n",&n)) sf_error("Lacking parameter!");
   if((n1!=n*(n+1)/2) || (n11!=n)) sf_error("Dimension mistake!");
   
    A=sf_floatalloc(n*n);
    d=sf_floatalloc(n);

    sf_floatread(A,n*(n+1)/2,in);
    sf_floatread(d,n,rhs);

    cholesky_dec(A,n);
    cholesky_solve(A, d, n);  /*d is both the input right hand side and the output solution */
    
    if(verb)
    {
	for(i=0;i<n;i++)
	sf_warning("d[%d]=%f",i+1,d[i]);
    }
    
   sf_putint(out,"n1",n);
   sf_putint(out,"n2",1);
   sf_floatwrite(d,n,out);

   free(A);
   free(d);
   exit(0);
}

