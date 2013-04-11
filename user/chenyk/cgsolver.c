/* Implement conjugate gradient algorithm to solve symmetric positive definite matrix linear equation */
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
/*^*/

#include "matoper.h"
#include "vecoper.h"
#include "cgsolver.h"
#ifndef EPS
#define EPS 0.0000000000000000001
#endif
static int j;

void cg_solve(	float **a /*matrix operator*/, 
	      	float *x /*initial vector and final output solution*/,
		float *b /*right hand side*/,
		int n /*matrix size*/,
		int iter /*iteration value*/)
/*< solver using conjugate gradient algorithm >*/
{
  float *y, *r1, *r2, *p;
  float beta, alpha;
  y=sf_floatalloc(n); /* y is a transition vector */
  r1=sf_floatalloc(n);
  r2=sf_floatalloc(n);
  p=sf_floatalloc(n);
  
  matmult_init(a);
  matmult_lop(false, false, n, n ,x, y);
  scalesum(b,y,r2,n,1,-1);
  
  for(j=0;j<iter;j++)
    {
      if(j==0)
	{scale(r2,p,n,1);}
      else
	{beta=dotmultsum(r2,r2,n)/(dotmultsum(r1,r1,n)+EPS);
	  scalesum(r2,p,p,n,1,beta);}  /* r1 correpsonds to r_(k-2) and r2 corresponds to r_(k-1) */
      
      matmult_lop(false,false,n,n,p,y);
      alpha=dotmultsum(r2,r2,n)/(dotmultsum(p,y,n)+EPS);
      scalesum(x,p,x,n,1,alpha);
      scale(r2,r1,n,1);
      scalesum(r2,y,r2,n,1,-alpha);	
    }
  free(y);
  free(r1);
  free(r2);
  free(p);
}




