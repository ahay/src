/* Implement classic iterative algorithm to solve linear equation 
(Jacobi iteration, Gauss-Seidel iteration, Successive Over Relaxation iteration) */
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
#include "classicsolver.h"

#ifndef EPS
#define EPS 0.0000000000000000001
#endif

static int i,j,k;

void jacobi_solve(	
		  float **a /*matrix operator*/, 
		  float *x /*initial vector and final output solution*/,
		  float *b /*right hand side*/,
		  int n /*matrix size*/,
		  int iter /*iteration value*/)
/*< solver using jacobi algorithm >*/
{
  float *y, *d;
  
  y=sf_floatalloc(n); /* y is a transition vector */
  d=sf_floatalloc(n);
  
  for(i=0;i<n;i++)
    {d[i]=a[i][i];a[i][i]=0;} /*compute A=B+D = (L+U)+D */
  scale(a[0],a[0],n*n,-1);
  
  matmult_init(a); /* Here A=(L+U) */
  for(k=0;k<iter;k++)
    {
      matmult_lop(false,false,n,n,x,y);
      
      vecdiv(y,d,x,n);
      vecdiv(b,d,y,n);
      scalesum(x,y,x,n,1,1);
    }
  
  scale(a[0],a[0],n*n,-1);
  for(i=0;i<n;i++) a[i][i]=d[i]; /*recovery for continuous computation*/
  free(y);
  free(d);
}

void gs_solve(	float **a /*matrix operator*/, 
	      	float *x /*initial vector and final output solution*/,
		float *b /*right hand side*/,
		int n /*matrix size*/,
		int iter /*iteration value*/)
/*< solver using gauss-seidel algorithm >*/
{
  float *d, *g, **aa, t;
  
  d=sf_floatalloc(n);
  g=sf_floatalloc(n);
  aa=sf_floatalloc2(n,n);
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)	
      {if(i!=j) aa[i][j]=-a[i][j]; else {aa[i][j]=0;d[i]=a[i][j];}} /*compute A=D-L-U */
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(i!=j)aa[i][j]=aa[i][j]/(d[i]+EPS); /*compute B=D\(A-D)*/
  
  vecdiv(b,d,g,n); /*compute g=D\b */
  for(k=0;k<iter;k++) 
    /* Because of the recursive property of Gauss-Seidel, 
   simple matrix or vector algebra operation should be 
   substituded by scaler operation. */
    {
      for(i=0;i<n;i++)
	{
	  t=0;
	  for(j=0;j<i;j++) t+=aa[i][j]*x[j];
	  for(j=i+1;j<n;j++) t+=aa[i][j]*x[j];
	  t+=g[i];
	  x[i]=t;
	}
    }
  
  free(aa[0]);
  free(aa);
  free(d);
}

void sor_solve(	float **a /*matrix operator*/, 
	      	float *x /*initial vector and final output solution*/,
		float *b /*right hand side*/,
		int n /*matrix size*/,
		int iter /*iteration value*/,
		float w /*relaxation coefficient*/)
/*< solver using succesive over relaxation algorithm>*/
{
  float *d, *g, **aa, t;
  
  aa=sf_floatalloc2(n,n);
  d=sf_floatalloc(n);
  g=sf_floatalloc(n);
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)	
      {if(i!=j) aa[i][j]=-a[i][j]; else {aa[i][j]=0;d[i]=a[i][j];}} /*compute A=D-L-U */
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(i!=j)aa[i][j]=aa[i][j]/(d[i]+EPS); /*compute B=D\(A-D)*/
  
  vecdiv(b,d,g,n); /*compute g=D\b */
  for(k=0;k<iter;k++) 
    /* Because of the recursive property of Gauss-Seidel and SOR, 
       simple matrix or vector algebra operation should be 
       substituded by scaler operation. */
    {
      for(i=0;i<n;i++)
	{
	  t=0;
	  t+=(1-w)*x[i];
	  for(j=0;j<i;j++) t+=w*aa[i][j]*x[j];
	  for(j=i+1;j<n;j++) t+=w*aa[i][j]*x[j];
	  t+=w*g[i];
	  x[i]=t;
	}
    }
  free(aa[0]);
  free(aa);
  free(d);
}


