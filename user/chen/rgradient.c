/* 3D rgradient by finite difference */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lphpoly.h"
#include "recursion.h"

static	float **c, ***g, **b, **b1, **b2;
static int order, n1, n2, nf;
static char mode;


void rgradient_init(char* type, int horder, int m1, int m2)
/*< initialize >*/
{
	float **p;
	int i;

	order = horder;
	nf = 2*order+1;
	c = sf_floatalloc2(nf, 2);
	p = lphpoly(order, order, type);
	for(i=0; i<2*nf; i++)	c[0][i] = p[0][i];
	free(p[0]);
	free(p);
	mode = type[0];

	n1 = m1;
	n2 = m2;

	b1 = sf_floatalloc2(n1, n2);
	b2 = sf_floatalloc2(n1, n2);

	g = sf_floatalloc3(n1*3, n2, nf);
	b = g[0];
	memset(g[0][0], 0, n1*n2*nf*sizeof(float));
#ifdef _OPENMP
    omp_init();
#endif
}

void rgradient_close()
/*< release the memory >*/
{
	free(c[0]);
	free(c);
	free(*b);
	free(b);
	free(*b1);
	free(b1);
	free(*b2);
	free(b2);
	free(g);
}

void rgradient(float **u1, float **u2)
/*< recursive computation of rgradient >*/
{
	int i1, i2, j1, j2, j3;
	float **pp;

#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1, i2, j1, j2, j3)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		b1[i2][i1] = 0.0;
		b2[i2][i1] = 0.0;
		if(i1-order<0 || i1+order>=n1) continue;	
		for(j1=-order; j1<=order; j1++)
		{
			b1[i2][i1] += c[1][j1+order]*u1[i2][i1+j1];
			b2[i2][i1] += c[0][j1+order]*u1[i2][i1+j1];
		}
	}


#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1, i2, j1, j2, j3)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		u2[i2][i1*3] = 0.0;
		u2[i2][i1*3+1] = 0.0;
		u2[i2][i1*3+2] = 0.0;
		if(i2-order<0 || i2+order>=n2) continue;	
		for(j2=-order; j2<=order; j2++)
		{
			u2[i2][i1*3] += c[0][j2+order]*b1[i2+j2][i1];
			u2[i2][i1*3+1] += c[1][j2+order]*b2[i2+j2][i1];
			u2[i2][i1*3+2] += c[0][j2+order]*b2[i2+j2][i1];
		}
	}

    pp = g[nf-1];
    for(i2=1; i2<nf; i2++) g[i2] = g[i2-1];
    g[0] = pp;
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1*3; i1++)
	g[0][i2][i1] = u2[i2][i1];


#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1, i2, j1, j2, j3)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		u2[i2][i1*3] = 0.0;
		u2[i2][i1*3+1] = 0.0;
		u2[i2][i1*3+2] = 0.0;
		for(j3=0; j3<nf; j3++)
		{
			u2[i2][i1*3] += c[0][j3] * g[j3][i2][i1*3];
			u2[i2][i1*3+1] += c[0][j3] * g[j3][i2][i1*3+1];
			u2[i2][i1*3+2] += c[1][j3] * g[j3][i2][i1*3+2];
		}
	}
}


