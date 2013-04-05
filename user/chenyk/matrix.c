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
/* multiplication: a->m*n, b->n*k */
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=k-1; j++)	{ 
	    u=i*k+j; c[u]=0.0;
	    for (l=0; l<=n-1; l++)
		c[u]+=a[i*n+l]*b[l*k+j];
	}
    return;
}

void add(float *a,float *b,int m,int n,float *c) 
/* add: a->m*n, b->m*n */
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]+b[i*n+j];
	}
    return;
}

void sub(float *a,float *b,int m,int n,float *c) 
/* subtract: a->m*n, b->m*n */
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]-b[i*n+j];
	}
    return;
}

void dotmul(float *a,float *b,int m,int n,float *c) 
/* dot multiplication: a->m*n, b->m*n */
{
    int i,j;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=n-1; j++)	{ 
		c[i*n+j]=a[i*n+j]*b[i*n+j];
	}
    return;
}


