/* Symmetric positive definite matrix equation solver using square root method (cholesky decomposition method) */
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
#include "cholesky.h"
static int i, k, j;

void cholesky_dec(	float *a/*input matrix and output triangle matrix*/,
			int n/*matrix size*/)
/*< Cholesky decomposition >*/
{
	for(k=1;k<=n;k++)
	{
		a[k*(k+1)/2-1]=sqrtf(a[k*(k+1)/2-1]);  /*A[i,j]=a[i*(i-1)/2+j-1], A[2,1]=a[1], where i and j (i>=j) are the index of matrix axis*/
		for(i=k+1;i<=n;i++)
			a[i*(i-1)/2+k-1]=a[i*(i-1)/2+k-1]/a[k*(k+1)/2-1];
		for(j=k+1;j<=n;j++)
			for(i=j;i<=n;i++)
				a[i*(i-1)/2+j-1]=a[i*(i-1)/2+j-1]-a[i*(i-1)/2+k-1]*a[j*(j-1)/2+k-1];
	}
}

void cholesky_solve(    float *l/*input triangle matrix*/, 
			float *d/*input right hand side and output solution*/, 
			int n/*matrix size*/)
/*< Solving equation >*/
{
	/* Solving down triangle matrix equation:Ly=d */
	d[0]=d[0]/l[0];
	for(i=2;i<=n;i++)
	{
		d[i-1]=d[i-1]/l[i*(i+1)/2-1];
		for(j=1;j<=i-1;j++)
			d[i-1]-=l[i*(i-1)/2+j-1]*d[j-1]/l[i*(i+1)/2-1];

	}
	
	/* Solving upper triangle matrix equation:L'x=y */
	d[n-1]=d[n-1]/l[n*(n+1)/2-1];
	for(i=n-1;i>=1;i--)
	{
		d[i-1]=d[i-1]/l[i*(i+1)/2-1];
		for(j=i+1;j<=n;j++)
			d[i-1]-=l[j*(j-1)/2+i-1]*d[j-1]/l[i*(i+1)/2-1];

	}
}


