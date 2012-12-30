/* Total Least Square eatimation */

/*
  Copyright (C) 2013 Zhonghuan Chen, UT Austin, Tsinghua University
  
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

struct tag_tls2{
	void *h;
	int par[8];
};

static void tls(float *out, float **in, int m1, int m2, int *par)
{
	int i1, i2, j1, j2, j3, n1, n2, rc1, rc2, rc3;
	float a, b, t1, t2, t3;
	n1 = par[0];
	n2 = par[1];
	rc1 = par[2];
	rc2 = par[3];
	rc3 = par[4];

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for(j3=0; j3<m2; j3++)
		for(j2=-rc2; j2<=rc2; j2++)
		for(j1=-rc1; j1<=rc1; j1++)
		{
			if(j1+i1<=0 || j1+i1>= n1) continue;
			if(j2+i2<=0 || j2+i2>= n2) continue;
			a = in[j3][(i2*n1+i1)*2];
			b = in[j3][(i2*n1+i1)*2+1];
			t1 += a*a;
			t2 += a*b;
			t3 += b*b;
		}
		a = t1 -t3;
		b = sqrt(a*a + 4*t2*t2);
		out[i2*n1+i1] = atan2(t2, (a-b)/2.0);
	}
}

static void ls(float *out, float **in, int m1, int m2, int *par)
{
	int i1, i2, j1, j2, j3, n1, n2, rc1, rc2, rc3;
	float a, b, t1, t2;
	n1 = par[0];
	n2 = par[1];
	rc1 = par[2];
	rc2 = par[3];
	rc3 = par[4];

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		t1 = 0.0;
		t2 = 0.0;
		for(j3=0; j3<m2; j3++)
		for(j2=-rc2; j2<=rc2; j2++)
		for(j1=-rc1; j1<=rc1; j1++)
		{
			if(j1+i1<=0 || j1+i1>= n1) continue;
			if(j2+i2<=0 || j2+i2>= n2) continue;
			a = in[j3][(i2*n1+i1)*2];
			b = in[j3][(i2*n1+i1)*2+1];
			t1 += a*a;
			t2 += a*b;
		}
		out[i2*n1+i1] = atan2(t2, t1);
	}
}

void* tls2_init(int m1, int m2, int *rect, bool istls)
/*< initialize >*/
{
	struct tag_tls2 *p;
	int n;

	p = sf_alloc(1,sizeof(struct tag_tls2));

	p->par[0] = m1;
	p->par[1] = m2;
	p->par[2] = rect[0];
	p->par[3] = rect[1];
	p->par[4] = rect[2];

	n = 2*rect[2]+1;
	if(istls) p->h =recursion_init(m1*m2*2, n, tls, p->par);
	else p->h =recursion_init(m1*m2*2, n, ls, p->par);
	return p;
}

void tls2_close(void *h)
/*< release memory >*/
{
	struct tag_tls2 *p = (struct tag_tls2*) h;
	recursion_close(p->h);
	free(h);
}

void tls2(void *h, float *u1)
/*< 2*n1*n2 ==> n1*n2 >*/
{
	struct tag_tls2 *p = (struct tag_tls2*) h;
	recursion(p->h, u1);
}

