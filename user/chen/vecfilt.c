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

static	float **buf, **b1, *b2, *d1, *d2;
static int n0, n1, n2, rect[3], n00;
static void *h3;


void vecfilt_mean(float *out, float **in, int m1, int m2, int *par)
{
	int i1, i2, j2, off;

	for(i2=0; i2<n1*n2; i2++)
	for(i1=0; i1<n00; i1++)
	{
		off = i2*(n00+1)+i1;
		out[off] = 0.0;
		for(j2=0; j2<m2; j2++) 
			out[off] += in[j2][off]*in[j2][i2*(n00+1)+n00];
		out[off] /= m2;
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

	n00 =m0*(m0+1)/2;

	if(rect[2]>0)
	h3 = recursion_init((n00+1)*m1*m2, 2*rect[2]+1, vecfilt_mean, NULL);

	buf = sf_floatalloc2(n00+1, m1*m2);
	b1 = sf_floatalloc2(m0, m0);
	b2 = sf_floatalloc(m0*m0);
	d1 = sf_floatalloc(n1);
	d2 = sf_floatalloc(n2);
}

void vecfilt_close()
/*< release the memory >*/
{
	free(buf[0]);
	free(buf);
	free(*b1);
	free(b1);
	free(b2);
	free(d1);
	free(d2);

	if(rect[2]>0)	recursion_close(h3);
}


void vecfilt(float *u1, float *u2, float *w)
/*< u1[n2][n1][3] u2[n2][n1][3] >*/
{
	int i0, i1, i2, inc=1, j1, lwork, info;
	float alpha=1.0, *q;

	if (rect[0]==0 && rect[1]==0 && rect[2] == 0)
	{
		for(i1=0; i1<n1*n2*n0; i1++) u2[i1] = u1[i1];
		return;
	}

	for(i2=0; i2<n1*n2; i2++)
	{
		for(i1=0; i1<n0*n0; i1++) b1[0][i1] = 0.0;
		ssyr_("U", &n0, &alpha, u1+i2*n0, &inc, b1[0], &n0); 
		for(i1=0, j1=0; i1<n0; i1++)
		for(i0=0; i0<=i1; i0++)
			buf[i2][j1++] = b1[i1][i0];
		if(w) buf[i2][n00] = w[i2];
		else buf[i2][n00] = 1.0;
	}

	if(rect[0]>0)
	{
		for(i2=0; i2<n2; i2++)
		for(i0=0; i0<n00; i0++)
		{
			for(i1=0; i1<n1; i1++)
			{
				d1[i1] = 0.0;
				for(j1=-rect[0]; j1<=rect[0]; j1++)
				{
					if(i1+j1<0) continue;
					if(i1+j1>=n1) break;
					d1[i1] += buf[i2*n1+i1+j1][i0] * buf[i2*n1+i1+j1][n00];
				}
			}
			for(i1=0; i1<n1; i1++)
				buf[i2*n1+i1][i0] =  d1[i1]/n1;
		}
	}
	if(rect[1]>0)
	{
		for(i1=0; i1<n1; i1++)
		for(i0=0; i0<n00; i0++)
		{
			for(i2=0; i2<n2; i2++)
			{
				d2[i2] = 0.0;
				for(j1=-rect[1]; j1<=rect[1]; j1++)
				{
					if(i2+j1<0) continue;
					if(i2+j1>=n2) break;
					d2[i2] += buf[(i2+j1)*n1+i1][i0] * buf[(i2+j1)*n1+i1][n00];
				}
			}
			for(i2=0; i2<n2; i2++)
				buf[i2*n1+i1][i0] =  d2[i2]/n2;
		}
	}

	if(rect[2]>0) recursion(h3, buf[0]);
	
	// pca
	lwork = n0*n0;

	for(i2=0; i2<n1*n2; i2++)
	{
		for(i1=0, j1=0; i1<n0; i1++)
		for(i0=0; i0<=i1; i0++)
		b1[i1][i0] = buf[i2][j1++];
		q = u2+i2*n0;
		ssyev_("V", "U", &n0, b1[0], &n0, q, b2, &lwork, &info);
		for(i0=0, j1=0; i0<n0; i0++)
		if(q[i0] > q[j1]) j1=i0;
		for(i0=0; i0<n0; i0++)
		q[i0] = b1[j1][i0];
	}
}



