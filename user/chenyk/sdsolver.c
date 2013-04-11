/* Implement steepest descent algorithm to solve symmetric positive definite matrix linear equation */
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
#include "sdsolver.h"
#ifndef EPS
#define EPS 0.0000000000000000001
#endif
static int j;

void sd_solve(	float **a /*matrix operator*/, 
	      	float *x /*initial vector and final output solution*/,
		float *b /*right hand side*/,
		int n /*matrix size*/,
		int iter /*iteration value*/)
/*< solver using steepest descent algorithm >*/
{
  float *y, *r;
  float alpha;
  y=sf_floatalloc(n); /* y is a transition vector */
  r=sf_floatalloc(n);
  
  matmult_init(a);
  matmult_lop(false, false, n, n ,x, y);
  scalesum(b,y,r,n,1,-1);
  
  for(j=0;j<iter;j++)
    {
      
      matmult_lop(false,false,n,n,r,y);
      alpha=dotmultsum(r,r,n)/(dotmultsum(r,y,n)+EPS);
      scalesum(x,r,x,n,1,alpha);
      matmult_lop(false,false,n,n,x,y);
      scalesum(b,y,r,n,1,-1);	
    }
  free(y);
  free(r);
}




