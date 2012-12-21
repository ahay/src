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
#include "dsp.h"
#include "rfir.h"

static	float **c, **b1, **b2;
static int order, n1, n2;
static char mode;
static void *h0, *h1;

void rgradient_init(char* type, int horder, int m1, int m2)
/*< initialize >*/
{
	float **p;
	int i, nf;

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

	h0 = rfir_init(nf, c[0], n1*n2);
	h1 = rfir_init(nf, c[1], n1*n2);

	b1 = sf_floatalloc2(n1, n2);
	b2 = sf_floatalloc2(n2, n1);

#ifdef _OPENMP
    omp_init();
#endif
}

void rgradient_close()
/*< release the memory >*/
{
	free(c[0]);
	free(c);
	free(*b1);
	free(*b2);
	free(b1);
	free(b2);
	rfir_close(h0);
	rfir_close(h1);
}


void rgradient(float *u1, float *u2)
/*< direct computation of rgradient >*/
{
	int i1, i2, off;

#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i2, off)
#endif
	for(i2=0; i2<n2; i2++)
	{
		off = i2*n1;
		firs(-order, order, c[1]+order, u1+off, 1, n1, u2+off*3, 3);
	}
#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1, off)
#endif
	for(i1=0; i1<n1; i1++)
	{
		firs(-order, order, c[1]+order, u1+i1, n1, n2, u2+i1*3+1, 3*n1);
	}
	rfir(h1, u1);
	for(i1=0; i1<n1*n2; i1++) u2[i1*3+2] = u1[i1];

	if(0)
	if(mode !='l')
	{
#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1,i2, off)
#endif
		for(i2=0; i2<n2; i2++) // d1
		{
			off = i2*n1;
			firs(-order, order, c[0]+order, u2+off*3+1, 3, n1, b1[i2], 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+1] = b1[i2][i1];
			firs(-order, order, c[0]+order, u2+off*3+2, 3, n1, b1[i2], 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+2] = b1[i2][i1];
		}
#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1,i2, off)
#endif
		for(i1=0; i1<n1; i1++) // d2
		{
			off = i1;
			firs(-order, order, c[0]+order, u2+off*3, 3*n1, n2, b2[i1], 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3] = b2[i1][i2];
			firs(-order, order, c[0]+order, u2+off*3+2, 3*n1, n2, b2[i1], 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3+2] = b2[i1][i2];
		}

		for(i1=0; i1<n1*n2; i1++) b1[0][i1] = u2[i1*3];
		rfir(h0, b1[0]);
		for(i1=0; i1<n1*n2; i1++) u2[i1*3] = b1[0][i1];
		for(i1=0; i1<n1*n2; i1++) b1[0][i1] = u2[i1*3+1];
		rfir(h0, b1[0]);
		for(i1=0; i1<n1*n2; i1++) u2[i1*3+1] = b1[0][i1];

	}
}

