/* Solver for complex linear equations  */
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

#include "csolver.h"
#include "solver.h"
#include "classicsolver.h"
#include "sdsolver.h"

static int i1,i2,n;
static float **aa;
static float *bb;
static float *xx;

void csolve_init(kiss_fft_cpx **a,kiss_fft_cpx *b,int N)
/*< initiate complex solver >*/
{
    n=N;
    aa=sf_floatalloc2(2*n,2*n);
    bb=sf_floatalloc(2*n);
    xx=sf_floatalloc(2*n);

    for(i2=0;i2<n;i2++)
	for(i1=0;i1<n;i1++)
	     {	aa[2*i2][2*i1]=a[i2][i1].r;
		aa[2*i2+1][2*i1+1]=a[i2][i1].r;
	      	aa[2*i2][2*i1+1]=-a[i2][i1].i;
		aa[2*i2+1][2*i1]=a[i2][i1].i;	}
    for(i1=0;i1<n;i1++)
	{	
	bb[2*i1]=b[i1].r;bb[2*i1+1]=b[i1].i;	
	}
    
 /*   for(i2=0;i2<2*n;i2++)
	{
	for(i1=0;i1<2*n;i1++)
	    printf("aa[%d][%d]=%4.2f",i2,i1,aa[i2][i1]); 
	printf("\n");
	}
    for(i1=0;i1<2*n;i1++)
	printf("bb[%d]=%4.2f",i1,bb[i1]);
    printf("\n");*/
}

void gaussel_csolve(kiss_fft_cpx *x, int niter)
/*< Gauss elimination complex solver >*/
{
    gaussel_init (2*n);
    gaussel_solve(aa, bb, xx);
    gaussel_close();

    for(i1=0;i1<n;i1++)
    {x[i1].r=xx[2*i1]; x[i1].i=xx[2*i1+1];}
}

void gs_csolve(kiss_fft_cpx *x, int niter)
/*< Gauss-Seidel iteration complex solver >*/
{
    gs_solve(aa, xx, bb, 2*n, niter);
    for(i1=0;i1<n;i1++)
    {x[i1].r=xx[2*i1]; x[i1].i=xx[2*i1+1];}
}

void jacobi_csolve(kiss_fft_cpx *x, int niter)
/*< Jacobi iteration iteration complex solver >*/
{
    jacobi_solve(aa, xx, bb, 2*n, niter);	
    for(i1=0;i1<n;i1++)
    {x[i1].r=xx[2*i1]; x[i1].i=xx[2*i1+1];}
}

void sor_csolve(kiss_fft_cpx *x, int niter, float w)
/*< SOR iteration iteration complex solver >*/
{
    sor_solve(aa, xx, bb, 2*n, niter, w);
    for(i1=0;i1<n;i1++)
    {x[i1].r=xx[2*i1]; x[i1].i=xx[2*i1+1];}
}

void sd_csolve(kiss_fft_cpx *x, int niter)
/*< Steepest descent complex solver >*/
{
    sd_solve(aa, xx, bb, 2*n, niter);
    for(i1=0;i1<n;i1++)
    {x[i1].r=xx[2*i1]; x[i1].i=xx[2*i1+1];}
}






