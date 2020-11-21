/* first coherence */

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

static int nw, n1, n2, lag1, lag2;
static float ***u0, ***u1, **v;

static void coh1_normalize(float *d, float **dn)
/* normalize a trace */
{
	int i1, iw;
	double t1, t2;

	for(iw=0, t2=0.0; iw<=2*nw; iw++)
		t2 += d[iw]*d[iw];
	t1 = sqrt(t2)+FLT_EPSILON;

	for(iw=0; iw<=2*nw; iw++)
	dn[nw][iw] = d[iw-nw]/t1;
	
	for(i1=nw+1; i1<n1-nw; i1++)
	{
		t1 = d[i1+nw];
		t2 += t1*t1; 
		t1 = d[i1-nw-1];
		t2 -= t1*t1; 
		t1 = sqrt(t2)+FLT_EPSILON;
		for(iw=-nw; iw<=nw; iw++)
		dn[i1][iw+nw] = d[iw+i1]/t1;
	}
	for(i1=0; i1<nw; i1++)
	for(iw=0; iw<=2*nw; iw++)
	{
		dn[i1][iw] = 0.0;
		dn[n1-i1-1][iw] = 0.0;
	}

}

static int dipscan(int n, int minl, int maxl, float *x, float **y, float *pv)
/* dip scan */
{
	int il, i1, max, nl;

	max = 0;
	nl = maxl-minl;
	for(il=0; il<=nl; il++)
	{
		pv[il] = 0.0;
		for(i1=0; i1<=2*n; i1++)
			pv[il] += x[i1]*y[il+minl][i1];
		if(pv[il] > pv[max]) max = il;
	}
	return max;
}

void coh1_init(int win, int m1, int m2, int l1, int l2)
/*< initialize >*/
{
	nw = win;
	n1 = m1;
	n2 = m2;
	lag1 = l1;
	lag2 = l2;
	u0 = sf_floatalloc3(2*nw+1, n1, n2);
	u1 = sf_floatalloc3(2*nw+1, n1, n2);
	v = sf_floatalloc2(2*(l1>l2? l1:l2)+1, n2);
	/* use two d memory for OPENMP */
}


void coh1_close()
/*< release memory >*/
{
	free(u0[0][0]);
	free(u1[0][0]);
	free(u0[0]);
	free(u1[0]);
	free(v[0]);
	free(u0);
	free(u1);
	free(v);
}


#define MY_MIN(a,b)  (a)<(b)?(a):(b)
#define MY_MAX(a,b)  (a)>(b)?(a):(b)
void coh1_2d(float **d)
/*< two-d coherence >*/
{
	int i1, i2, k;

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,8)   \
	private(i2)       
#endif
	for(i2=0; i2<n2; i2++)
		coh1_normalize(d[i2], u1[i2]);

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,1)   \
	private(i1,i2,k)       
#endif
	for(i2=1; i2<n2; i2++)
	{
		for(i1=0; i1<n1; i1++)
		{
			k = dipscan(nw, MY_MAX(-i1, -lag1), MY_MIN(n1-i1-1, lag1),
				u1[i2-1][i1], u1[i2]+i1, v[i2]);
			d[i2-1][i1] = sqrt(v[i2][k]);
		}
	}
}


void coh1_3d_init2d(float **d)
/*< initialize for each 2D profile  >*/
{
	int i2;
#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,8)   \
	private(i2)       
#endif
	for(i2=0; i2<n2; i2++)
		coh1_normalize(d[i2], u0[i2]);
}


void coh1_3d(float **d)
/*< recursive first coherence >*/
{
	int i1, i2, k1, k2;
	float ***p, t1;

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,8)   \
	private(i2)       
#endif
	for(i2=0; i2<n2; i2++)
		coh1_normalize(d[i2], u1[i2]);

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,1)   \
	private(i1,i2, t1, k1, k2)       
#endif
	for(i2=1; i2<n2; i2++)
	{
		for(i1=0; i1<n1; i1++)
		{
			k1 = dipscan(nw, MY_MAX(-i1, -lag1), MY_MIN(n1-i1-1, lag1),
				u1[i2-1][i1], u1[i2]+i1, v[i2]);
			t1 = v[i2][k1];
			k2 = dipscan(nw, MY_MAX(-i1, -lag2), MY_MIN(n1-i1-1, lag2),
				u1[i2-1][i1], u0[i2-1]+i1, v[i2]);
			d[i2-1][i1] = sqrt(t1*v[i2][k2]);
		}
	}
	p = u0;
	u0 = u1;
	u1 = p;
}



