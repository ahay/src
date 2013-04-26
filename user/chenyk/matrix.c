/* Matrix algebra operation*/
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>
#include "matrix.h"

void mul(float *a,float *b,int m,int n,int k,float *c) 
/*< multiplication: a->m*n, b->n*k >*/
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=k-1; j++)	{ 
	    u=i*k+j; c[u]=0.0;
	    for (l=0; l<=n-1; l++)
		c[u]+=a[i*n+l]*b[l*k+j];
	}
}

void add(float *a,float *b,int m,int n,float *c) 
/*< add: a->m*n, b->m*n >*/
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]+b[i*n+j];
	}
}

void sub(float *a,float *b,int m,int n,float *c) 
/*< subtract: a->m*n, b->m*n >*/
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]-b[i*n+j];
	}
}

void dotmul(float *a,float *b,int m,int n,float *c) 
/*< dot multiplication: a->m*n, b->m*n >*/
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]*b[i*n+j];
	}
}

void tran(float *a, int m, int n, float *b)
/*< matrix transpose: a->m*n, b->n*m >*/
{
    int i, j;
    for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		b[i*m+j]=a[j*n+i];
}

void zero(float *a, int m)
/*< zero a specific part of an array, like initialization >*/
{
    int i;
    for(i=0;i<m;i++)
	{ a[i]=0; }   
}

void cmatmul(kiss_fft_cpx *a,kiss_fft_cpx *b,int m,int n,int k,kiss_fft_cpx *c) 
/*< complex matrix multiplication: a->m*n, b->n*k >*/
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=k-1; j++)	
	{ 
	    u=i*k+j; c[u].r=0.0; c[u].i=0.0;
	    for (l=0; l<=n-1; l++)
		c[u]=sf_cadd(c[u],sf_cmul(a[i*n+l],b[l*k+j]));
	}
}

void cmattran(kiss_fft_cpx *a, int m, int n, kiss_fft_cpx *b)
/*< complex matrix transpose: a->m*n, b->n*k >*/
{
    int i, j;
    for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		b[i*m+j]=a[j*n+i];	

}

void czero(kiss_fft_cpx *a, int m )
/*< zero a specific part of a complex array, like initialization >*/
{
    int i;
    for(i=0;i<m;i++)
	{a[i].r=0;a[i].i=0;}   
}



