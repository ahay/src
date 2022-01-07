/* first coherence with estimated dip */

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

static int nw, n1, n2;
static float ***u0, ***u1; /* **v; */

static void coh1dip_normalize(float *d, float **dn)
// normalize a trace
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


static float coh1dip(int n, float *x, float *y)
{
	int i1;
	float p;

	p = 0.0;
	for(i1=0; i1<=2*n; i1++)
		p += x[i1]*y[i1];

	return p;
}


void coh1dip_init(int win, int m1, int m2)
/*< initialize >*/
{
	nw = win;
	n1 = m1;
	n2 = m2;
	u0 = sf_floatalloc3(2*nw+1, n1, n2);
	u1 = sf_floatalloc3(2*nw+1, n1, n2);
	// use two d memory for OPENMP
}


void coh1dip_close()
/*< release memory >*/
{
	free(u0[0][0]);
	free(u1[0][0]);
	free(u0[0]);
	free(u1[0]);
	free(u0);
	free(u1);
}


#define MY_MIN(a,b)  (a)<(b)?(a):(b)
#define MY_MAX(a,b)  (a)>(b)?(a):(b)
void coh1dip_2d(float **d, float **p)
/*< two-d coherence >*/
{
	int i1, i2, k;

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,8)   \
	private(i2)       
#endif
	for(i2=0; i2<n2; i2++)
		coh1dip_normalize(d[i2], u1[i2]);

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,1)   \
	private(i1,i2,k)       
#endif
	for(i2=1; i2<n2; i2++)
	{
		for(i1=0; i1<n1; i1++)
		{
			k = i1 + p[i2-1][i1]+0.5;
			k = MY_MAX(k, 0);
			k = MY_MIN(k, n1-1);
			d[i2-1][i1] = sqrt(coh1dip(nw, u1[i2-1][i1], u1[i2][k]));
		}
	}
}


void coh1dip_3d_init2d(float **d)
/*< initialize for each 2D profile  >*/
{
	int i2;
#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,8)   \
	private(i2)       
#endif
	for(i2=0; i2<n2; i2++)
		coh1dip_normalize(d[i2], u0[i2]);
}


void coh1dip_3d(float **d, float **p1, float **p2)
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
		coh1dip_normalize(d[i2], u1[i2]);

#ifdef _OPENMP
#pragma omp parallel for     \
	schedule(dynamic,1)   \
	private(i1,i2, t1,k1, k2)       
#endif
	for(i2=1; i2<n2; i2++)
	{
		for(i1=0; i1<n1; i1++)
		{
			k1 = i1 + p1[i2-1][i1]+0.5;
			k1 = MY_MAX(k1, 0);
			k1 = MY_MIN(k1, n1-1);
			k2 = i1 + p2[i2-1][i1]+0.5;
			k2 = MY_MAX(k2, 0);
			k2 = MY_MIN(k2, n1-1);
			t1 = coh1dip(nw, u1[i2-1][i1], u1[i2][k1]);
			t1 *= coh1dip(nw, u0[i2-1][i1], u1[i2-1][k1]);
			d[i2-1][i1] = sqrt(t1);
		}
	}
	p = u0;
	u0 = u1;
	u1 = p;
}



