/* 3D rgradient by finite difference */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

static	float **c, *buf;
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

	buf = sf_floatalloc(n1*n2);

#ifdef _OPENMP
    omp_init();
#endif
}

void rgradient_close()
/*< release the memory >*/
{
	free(c[0]);
	free(c);
	free(buf);
	rfir_close(h0);
	rfir_close(h1);
}


void rgradient(float *u1, float *u2)
/*< direct computation of rgradient >*/
{
	int i1, i2, off;

	for(i2=0; i2<n2; i2++)
	{
		off = i2*n1;
		firs(order, order, c[1], u1+off, 1, n1, u2+off*3, 3);
	}
	for(i1=0; i1<n1; i1++)
	{
		firs(order, order, c[1], u1+i1, n1, n2, u2+i1*3+1, 3*n1);
	}
	rfir(h1, u1);
	for(i1=0; i1<n1*n2; i1++) u2[i1*3+2] = u1[i1];

	if(mode !='l')
	{
		for(i2=0; i2<n2; i2++) // d1
		{
			off = i2*n1;
			firs(order, order, c[0], u2+off*3+1, 3, n1, buf, 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+1] = buf[i1];
			firs(order, order, c[0], u2+off*3+2, 3, n1, buf, 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+2] = buf[i1];
		}
		for(i1=0; i1<n1; i1++) // d2
		{
			off = i1;
			firs(order, order, c[0], u2+off*3, 3*n1, n2, buf, 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3] = buf[i2];
			firs(order, order, c[0], u2+off*3+2, 3*n1, n2, buf, 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3+2] = buf[i2];
		}

		for(i1=0; i1<n1*n2; i1++) buf[i1] = u2[i1*3];
		rfir(h0, buf);
		for(i1=0; i1<n1*n2; i1++) u2[i1*3] = buf[i1];
		for(i1=0; i1<n1*n2; i1++) buf[i1] = u2[i1*3+1];
		rfir(h0, buf);
		for(i1=0; i1<n1*n2; i1++) u2[i1*3+1] = buf[i1];

	}
}

