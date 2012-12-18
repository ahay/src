/* recursive vector filter */

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
#include "recursion.h"
#include "_lapack.h"

static	float **buf, *buf2, *wgt;
static int n0, n1, n2, rect[3];
static void *h3;

void vecfilt_mean(float *out, float **in, int m1, int m2)
{
	int i1, i2;
	for(i1=0; i1<m1; i1++)
	{
		out[i1] = 0.0;
		if(wgt) 
			for(i2=0; i2<m2; i2++) out[i1] += in[i2][i1]*wgt[i2];
		else for(i2=0; i2<m2; i2++) out[i1] += in[i2][i1];
		out[i1] /= m2;
	}
}

void vecfilt_init(int m0, int m1, int m2, int *rc)
/*< initialize >*/
{
	n0 = m0;
	n1 = m1;
	n2 = m2;
	rect[0] = rc[0];
	rect[1] = rc[1];
	rect[2] = rc[2];

	if(rect[2]>0)
	h3 = recursion_init(m0*m0*m1*m2, 2*rect[2]+1, vecfilt_mean);

	buf = sf_floatalloc2(m0*m0, m1*m2);
	buf2 = sf_floatalloc((m1>m2?m1:m2)+(rect[0]+rect[1])*2);
#ifdef _OPENMP
    omp_init();
#endif
}

void vecfilt_close()
/*< release the memory >*/
{
	free(buf[0]);
	free(buf);
	free(buf2);
	
	if(rect[2]>0)	recursion_close(h3);
}


void vecfilt(float *u1, float *u2, float *w)
/*< u1[n2][n1][3] u2[n2][n1][3] >*/
{
	int i0, i1, i2, n00, inc=1, j1, lwork, info;
	float alpha=1.0;

	n00 = n0*n0;
#ifdef _OPENMP
#pragma omp parallel for       \
    schedule(dynamic,5)         \
    private(i1)
#endif
	for(i1=0; i1<n1*n2; i1++)
		ssyr_("U", &n0, &alpha, u1+i1*n0, &inc, buf[i1], &n0); 

	wgt = w;

	if(rect[2]>0) recursion(h3, buf[0]);

	if(rect[0]>1)
	for(i2=0; i2<n2; i2++)
	for(i0=0; i0<n00; i0++)
	{
		for(i1=0; i1<n1; i1++)
		{
			buf2[i1] = 0.0;
			for(j1=i1-rect[0]; j1<i1+rect[0]; j1++)
			{
				if(j1<0) continue;
				if(j1>=n1) break;
				if(w) buf2[i1] += buf[i2*n1+j1][i0]*w[i2*n1+j1];
				else buf2[i1] += buf[i2*n1+j1][i0];
			}
		}
		for(i1=0; i1<n1; i1++)
			buf[i2*n1+i1][i0] =  buf2[i1]/n1;
	}

	for(i1=0; i1<n1; i1++)
	for(i0=0; i0<n00; i0++)
	{
		for(i2=0; i2<n2; i2++)
		{
			buf2[i2] = 0.0;
			for(j1=i2-rect[1]; j1<i2+rect[1]; j1++)
			{
				if(j1<0) continue;
				if(j1>=n2) break;
				if(w) buf2[i2] += buf[i2*n1+j1][i0]*w[i2*n1+j1];
				else buf2[i2] += buf[i2*n1+j1][i0];
			}
		}
		for(i2=0; i2<n2; i2++)
			buf[i2*n1+i1][i0] =  buf2[i2]/n2;
	}
	
	// pca
	lwork = n1+n2;

	for(i1=0; i1<n1*n2; i1++)
	{
		ssyev_("V", "U", &n0, buf[i1], &n0, u2+i1*n0, buf2, &lwork, &info);
		for(i0=0, i2=0; i0<n0; i0++)
		if(buf2[i0]>buf2[i2]) i2=i0;
		for(i0=0; i0<n0; i0++)
		u2[i1*n0+i0] = buf[i1][i2*n0+i0];
	}
}

