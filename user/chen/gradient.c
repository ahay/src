/* 3D gradient by finite difference */

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

static	float **c;
static int order, rect[3];
static char mode;


void gradient_init(char* type, int horder, int *rc)
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

	rect[0] = rc[0];
	rect[1] = rc[1];
	rect[2] = rc[2];

#ifdef _OPENMP
    omp_init();
#endif
}

void gradient_close()
/*< release the memory >*/
{
	free(c[0]);
	free(c);
}


void gradient(float *u1, float *u2, int n1, int n2, int n3)
/*< direct computation of gradient >*/
{
	int i1, i2, i3, off;
	float *buf;

	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	{
		off = (i3*n2+i2)*n1;
		firs(order, order, c[1], u1+off, 1, n1, u2+off*3, 3);
	}
	for(i3=0; i3<n3; i3++)
	for(i1=0; i1<n1; i1++)
	{
		off = i3*n2*n1+i1;
		firs(order, order, c[1], u1+off, n1, n2, u2+off*3+1, 3*n1);
	}
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		off = i2*n1+i1;
		firs(order, order, c[1], u1+off, n1*n2, n3, u2+off*3+2, 3*n2*n1);
	}
	if(mode !='l')
	{
		off = n1>n2 ? n1 : n2;
		off = off>n3 ? off : n3;
		buf = sf_floatalloc(off);

		for(i3=0; i3<n3; i3++)
		for(i2=0; i2<n2; i2++)
		{
			off = (i3*n2+i2)*n1;
			firs(order, order, c[0], u2+off*3+1, 3, n1, buf, 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+1] = buf[i1];
			firs(order, order, c[0], u2+off*3+2, 3, n1, buf, 1);
			for(i1=0; i1<n1; i1++) u2[(off+i1)*3+2] = buf[i1];
		}
		for(i3=0; i3<n3; i3++)
		for(i1=0; i1<n1; i1++)
		{
			off = i3*n2*n1+i1;
			firs(order, order, c[0], u2+off*3, 3*n1, n2, buf, 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3] = buf[i2];
			firs(order, order, c[0], u2+off*3+2, 3*n1, n2, buf, 1);
			for(i2=0; i2<n2; i2++) u2[(off+i2*n1)*3+2] = buf[i2];
		}
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
		{
			off = i2*n1+i1;
			firs(order, order, c[0], u2+off*3, 3*n1*n2, n3, buf, 1);
			for(i3=0; i3<n3; i3++) u2[(off+i3*n1*n2)*3] = buf[i3];
			firs(order, order, c[0], u2+off*3+1, 3*n1*n2, n3, buf, 1);
			for(i3=0; i3<n3; i3++) u2[(off+i3*n1*n2)*3+1] = buf[i3];
		}

		free(buf);
	}
}

