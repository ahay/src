/* Tridiagonal matrix solver using chasing method */
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
#include<rsf.h>
#include "tridsolver.h"
static int i, n;
static float *a, *b, *c;

void trid_init(int N/* matrix size */,
		       float *A/*input matrix*/)
/*< initialize >*/
/*|ac   |    */		/*alpha		    *//* 1 gamma 	*/
/*|bac  |    */		/* beta alpha	    *//*     1 gamma   	*/
/*| bac |x=d */  /*->*/ /*       beta alpha *//*        1 gamma	*/
/*|  ...|    */ 	/*	       beta *//*	    1	*/
/*|   ba|    */		/*		    *//*	  	*/
{
    n=N;
    a=sf_floatalloc(n);
    b=sf_floatalloc(n);
    c=sf_floatalloc(n);
    for(i=0;i<n;i++)
	{a[i]=A[i*n+i];}
    for(i=1;i<n;i++)
	{b[i]=A[i*n+i-1];}
    for(i=0;i<n-1;i++)
	{c[i]=A[i*n+i+1];}	
}

void trid_solve( float *d /* in - right-hand side */, 
			   float *x /* out - solution */)
/*< invert the matrix >*/
{
    float *alpha, *beta, *gamma, *y;
    alpha=sf_floatalloc(n);
    beta=sf_floatalloc(n);
    gamma=sf_floatalloc(n);
    y=sf_floatalloc(n);

    /* LU decomposition */
    gamma[0]=c[0]/(a[0]+0.000000000000001);
    alpha[0]=a[0];
    for(i=1;i<n-1;i++)
	{alpha[i]=a[i]-b[i]*gamma[i-1];
    	gamma[i]=c[i]/(alpha[i]+0.000000000000001);}
    alpha[n-1]=a[n-1]-b[n-1]*gamma[n-2];
    for(i=1;i<n;i++)beta[i]=b[i];

    /* Solving Ly=d */    
    y[0]=d[0]/(alpha[0]+0.000000000000001);
    for(i=1;i<n;i++)
	y[i]=(d[i]-b[i]*y[i-1])/(alpha[i]+0.000000000000001);

    /* Solving Ux=y */ 
    x[n-1]=y[n-1];
    for(i=n-2;i>=0;i--)
	x[i]=y[i]-gamma[i]*x[i+1];

}

void trid_close(void)
/*< free allocated memory >*/
{ 
    free(a);
    free(b);
    free(c);
}

