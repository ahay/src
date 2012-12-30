/* recursive vector filter */

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
#include "recursion.h"
#include "_lapack.h"

static	float **buf, *buf2, **b1, **b2, *wgt;
static int n0, n1, n2, rect[3];
static void *h3;


void vecfilt_mean(float *out, float **in, int m1, int m2, int *par)
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
	h3 = recursion_init(m0*m0*m1*m2, 2*rect[2]+1, vecfilt_mean, NULL);

	buf = sf_floatalloc2(m0*m0, m1*m2);
	b1 = sf_floatalloc2(m1, m2);
	b2 = sf_floatalloc2(m2, m1);
	buf2 = sf_floatalloc((m1>m2?m1:m2)+(rect[0]+rect[1])*2);
}

void vecfilt_close()
/*< release the memory >*/
{
	free(buf[0]);
	free(buf);
	free(buf2);
	free(*b1);
	free(*b2);
	free(b1);
	free(b2);

	if(rect[2]>0)	recursion_close(h3);
}


void vecfilt(float *u1, float *u2, float *w)
/*< u1[n2][n1][3] u2[n2][n1][3] >*/
{
	int i0, i1, i2, n00, inc=1, j1, lwork, info, off;
	float alpha=1.0, *p, *q;

	if (rect[0]==0 && rect[1]==0 && rect[2] == 0)
	{
		for(i1=0; i1<n1*n2*n0; i1++) u2[i1] = u1[i1];
		return;
	}

	n00 = n0*n0;
	for(i1=0; i1<n1*n2; i1++)
		ssyr_("U", &n0, &alpha, u1+i1*n0, &inc, buf[i1], &n0); 

	wgt = w;

	if(rect[2]>0) recursion(h3, buf[0]);

	if(rect[0]>0)
	{
		for(i2=0; i2<n2; i2++)
		for(i0=0; i0<n00; i0++)
		{
			if(i0/n0 < i0%n0) continue;
			for(i1=0; i1<n1; i1++)
			{
				b1[i2][i1] = 0.0;
				for(j1=i1-rect[0]; j1<=i1+rect[0]; j1++)
				{
					if(j1<0) continue;
					if(j1>=n1) break;
					if(w) b1[i2][i1] += buf[i2*n1+j1][i0]*w[i2*n1+j1];
					else b1[i2][i1] += buf[i2*n1+j1][i0];
				}
			}
			for(i1=0; i1<n1; i1++)
				buf[i2*n1+i1][i0] =  b1[i2][i1]/n1;
		}
	}
	if(rect[1]>0)
	{
		for(i1=0; i1<n1; i1++)
		for(i0=0; i0<n00; i0++)
		{
			if(i0/n0 < i0%n0) continue;
			for(i2=0; i2<n2; i2++)
			{
				b2[i1][i2] = 0.0;
				for(j1=i2-rect[1]; j1<=i2+rect[1]; j1++)
				{
					if(j1<0) continue;
					if(j1>=n2) break;
					if(w) b2[i1][i2] += buf[j1*n1+i1][i0]*w[j1*n1+i1];
					else b2[i1][i2] += buf[j1*n1+i1][i0];
				}
			}
			for(i2=0; i2<n2; i2++)
				buf[i2*n1+i1][i0] =  b2[i1][i2]/n2;
		}
	}
	
	// pca
	lwork = n1;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		p = buf[i2*n1+i1];
		q = u2+(i2*n1+i1)*n0;
		ssyev_("V", "U", &n0, p, &n0, q, b1[i2], &lwork, &info);
		for(i0=0, j1=0; i0<n0; i0++)
		if(q[i0] > q[j1]) j1=i0;
		for(i0=0; i0<n0; i0++)
			q[i0] = p[j1*n0+i0];
	}
}



