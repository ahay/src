/* Azimuth CURVature */

/*
  Copyright (C) 2013 Zhonghuan Chen
  
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
#include "lphpoly.h"
#include "tls2.h"
#include "recursion.h"
#include "dsp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static int order, n1, n2, n3, rect[3];
static float *c1, *c2;
static int ind, na;
static float ***dd, *num, *den;
static float **u, **ut, **ux, **uy, **u2tx, **u2ty, **u2xx, **u2yy, **u2xy;
static void **h, *h1;


void acurv_der(float *d, float **buf, int m1, int m2, int *par)
{
	int i1, i2, j1;
#ifdef _OPENMP
#pragma omp parallel for                    \
    schedule(dynamic,10)         
#endif
	for(i2=0; i2<n2; i2++)
		firs(-order, order, c1+order, buf[order]+i2*n1, 1, n1, ut[i2], 1);
#ifdef _OPENMP
#pragma omp parallel for                    \
    schedule(dynamic,10)         
#endif
	for(i1=0; i1<n1; i1++)
	{
		firs(-order, order, c1+order, buf[order]+i1, n1, n2, ux[0]+i1, n1);
		firs(-order, order, c2+order, buf[order]+i1, n1, n2, u2xx[0]+i1, n1);
		firs(-order, order, c1+order, ut[0]+i1, n1, n2, u2tx[0]+i1, n1);
	}

	for(i1=0; i1<m1; i1++)
	{
		uy[0][i1] = 0.0;
		u2yy[0][i1] = 0.0;
		for(j1=0; j1<m2; j1++)
		{
			uy[0][i1] += buf[j1][i1]*c1[m2-j1];
			u2yy[0][i1] += buf[j1][i1]*c2[m2-j1];
		}
	}
#ifdef _OPENMP
#pragma omp parallel for                    \
    schedule(dynamic,10)         
#endif
	for(i2=0; i2<n2; i2++)
		firs(-order, order, c1+order, uy[i2], 1, n1, u2ty[i2], 1);
}

void acurv_init(
	int in1, int in2, int in3,
	int horder, char*interp,
	int* rc, int nazimuth)
/*< initialiaze >*/
{
	int nf, i;
	float **p;

	order = horder;
	n1 = in1;
	n2 = in2;
	n3 = in3;
	rect[0] = rc[0];
	rect[1] = rc[1];
	rect[2] = rc[2];
	ind = 0;
	na = nazimuth;

	nf = 2*order+1;
    c1 = sf_floatalloc(nf);
    c2 = sf_floatalloc(nf);
    p = lphpoly(order, order, interp);
    for(i=0; i<nf; i++)
	{
		c1[i] = p[0][i];	
		c2[i] = p[1][i];	
	}
	free(p[0]); free(p);

	dd = sf_floatalloc3(n1, n2, 9);
	u = dd[0];
	ut = dd[1];
	ux = dd[2];
	uy = dd[3];
	u2tx = dd[4];
	u2ty = dd[5];
	u2xx = dd[6];
	u2yy = dd[7];
	u2xy = dd[8];
	num = sf_floatalloc(n1*n2);
	den = sf_floatalloc(n1*n2);

	h1 = recursion_init(n1*n2, nf, acurv_der, NULL);

	h = sf_alloc(na, sizeof(void*));
    for(i=0; i<na; i++)
		h[i] = tls2_init(n1, n2, rect, false, false);
}

void acurv_close()
/*< release memory >*/
{
	int i;

	free(c1);
	free(c2);
	free(num);
	free(den);
	
	free(dd[0][0]);
	free(dd[0]);
	free(dd);

	recursion_close(h1);
    for(i=0; i<na; i++)
		tls2_close(h[i]);
}

int acurv(float **d)
/*< 3d recursive curvature >*/
{
//  time sequences for recursive operators:
//          0   od   od+rc        n3  n3+od  n3+od+rc
//  i3      |----|-----|-----------|----|-----|
//  read    |----------------------|
//  ut           |----------------------|
//  ux           |----------------------|
//  u2xt         |----------------------|
//  u2x2         |----------------------|
//  uy           |----------------------|
//  u2y2         |----------------------|
//  u2yt         |----------------------|
//  ur           |----------------------|
//  u2r2         |----------------------|
//  u2rt         |----------------------|
//  den          |----------------------|
//  num          |----------------------|
//  write              |----------------------|
	int n12, i12, i3;
	float theta, c, s, ur, u2tr, u2rr;

	n12 = n1*n2;
	for(i12=0; i12<n12; i12++)	u[0][i12] = d[0][i12];
	if(ind==0)
		for(i3=0; i3<order; i3++) recursion_push(h1, u[0]);
	else if(ind>=1 && ind <order) recursion_push(h1, u[0]);
	else if(ind >= order && ind < n3+order)
	{
		recursion(h1, u[0]);
		for(i3=0; i3<na; i3++)
		{
			theta = SF_PI/(na-1)*i3;
			c = cos(theta);
			s = sin(theta);
			for(i12=0; i12<n12; i12++)
			{
				ur = ux[0][i12]*c + uy[0][i12]*s;
				u2tr = u2tx[0][i12]*c + u2ty[0][i12]*s;
				u2rr = u2xx[0][i12]*c*c + u2yy[0][i12]*s*s
					 + u2xy[0][i12]*c*s;
				den[i12] = ut[0][i12]*ut[0][i12];
				num[i12] = (u2tr*ur-u2rr*ut[0][i12])*den[i12];
				den[i12] = pow(ur*ur+den[i12], 1.5)*fabs(ut[0][i12]);
			}
			for(i12=0; i12<n12; i12++)
				d[i3][i12] = num[i12];
			tls2(h[i3], d[i3], den);

/*			for(i12=0; i12<n12; i12++)
			{
				d[0][i12] = u[0][i12];
				d[1][i12] = ut[0][i12];
				d[2][i12] = ux[0][i12];
				d[3][i12] = uy[0][i12];
				d[4][i12] = u2tx[0][i12];
				d[5][i12] = u2ty[0][i12];
				d[6][i12] = u2xx[0][i12];
				d[7][i12] = u2yy[0][i12];
				d[8][i12] = num[i12];
				d[9][i12] = den[i12];
			}
*/		} 
	}else if(ind >= n3+order)
	for(i3=0; i3<na; i3++)
	{
		for(i12=0; i12<n12; i12++)
			d[i3][i12] = num[i12];
		tls2(h[i3], d[i3], den);
	}

	ind++;
	return ind - 1;
}



