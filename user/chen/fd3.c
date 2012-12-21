/* 3D finite difference */

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

static	double c0, c11, c12, c21, c22, c31, c32;

void fd3_init(float d1, float d2, float d3)
/*< initialize >*/
{
	double t;

	c11 = 0.0;	c12 = 0.0;
	c21 = 0.0;	c22 = 0.0;
	c31 = 0.0;	c32 = 0.0;

	if(d1 != 0.0)
	{
		t = d1;
		t = 1.0/(t*t);
		c11 = 4.0*t/3.0;
		c12=  -t/12.0;
	}

	if(d2 != 0.0)
	{
		t = d2;
		t = 1.0/(t*t);
		c21 = 4.0*t/3.0;
		c22=  -t/12.0;
	}
	if(d3 != 0.0)
	{
		t = d3;
		t = 1.0/(t*t);
		c31 = 4.0*t/3.0;
		c32=  -t/12.0;
	}
	c0  = -2.0 * (c11 + c12 + c21 + c22 +c31 + c32);
#ifdef _OPENMP
    omp_init();
#endif
}

void fd3_laplacian(int n1, int n2, int n3, float *uin, float *uout)
/*< Laplacian operator, 4th-order finite-difference >*/
{
	int i1, i2, i3, i, n12;
	double u;

	n12 = n1*n2;

#ifdef _OPENMP
#pragma omp parallel for                    \
	schedule(dynamic,n1)         \
	private(i1,i2,i3,i,u)
/*	shared(n1,n2,n3,n12,i,u,uin,uout,c0,c11,c12,c21,c22,c31,c32) */
#endif
	for (i3=0; i3 < n3; i3++) 
	for (i2=0; i2 < n2; i2++) 
	for (i1=0; i1 < n1; i1++) 
	{
		i = i3*n12+i2*n1+i1;
		u = c0*uin[i];
		if(n1 != 1 )
		{
			if(i1 >= 1) u += c11*uin[i-1];
			if(i1 >= 2) u += c12*uin[i-2];
			if(i1 < n1-1) u += c11*uin[i+1];
			if(i1 < n1-2) u += c12*uin[i+2];
		}
		if(n2 != 1 )
		{
			if(i2 >= 1) u += c21*uin[i-n1];
			if(i2 >= 2) u += c22*uin[i-2*n1];
			if(i2 < n2-1) u += c21*uin[i+n1];
			if(i2 < n2-2) u += c22*uin[i+2*n1];
		}
		if(n3 != 1 )
		{
			if(i3 >= 1) u += c31*uin[i-n12];
			if(i3 >= 2) u += c32*uin[i-2*n12];
			if(i3 < n1-1) u += c31*uin[i+n12];
			if(i3 < n1-2) u += c32*uin[i+2*n12];
		}
		uout[i] = u;
	}
}



