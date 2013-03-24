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
#ifdef _OPENMP
#include <omp.h>
#endif

struct tag_tls2{
	void *h;
	int n1;
	int n2;
	int rc1;
	int rc2;
	int rc3;
	float **u1;
	float **u2;
	int n0;
	bool ang;
};

static void ls_mean3(float *out, float **in, int m1, int m2, int *par)
{
	int i1, i2;

#ifdef _OPENMP
#pragma omp parallel for                    \
	schedule(dynamic,10)         \
	private(i1,i2)
#endif
	for(i1=0; i1<m1; i1++)
	for(out[i1]=0.0, i2=0; i2<m2; i2++)
	{
		out[i1] += in[i2][i1];
	}
}

void* tls2_init(int m1, int m2, int *rect, bool istls, bool ang)
/*< initialize >*/
{
	struct tag_tls2 *p;
	int n;

	p = sf_alloc(1,sizeof(struct tag_tls2));

	p->n1 = m1;
	p->n2 = m2;
	p->rc1 = rect[0];
	p->rc2 = rect[1];
	p->rc3 = rect[2];
	p->ang = ang;
	if(istls) p->n0 = 3;
	else p->n0 = 2;

	n = 2*rect[2]+1;
	p->h = recursion_init(m1*m2*p->n0, n, ls_mean3, NULL);
	p->u1 = sf_floatalloc2(m1*p->n0, m2);
	p->u2 = sf_floatalloc2(m1*p->n0, m2);

#ifdef _OPENMP
    omp_init();
#endif
	return p;
}

void tls2_close(void *h)
/*< release memory >*/
{
	struct tag_tls2 *p = (struct tag_tls2*) h;
	recursion_close(p->h);
	free(p->u1[0]);
	free(p->u1);
	free(p->u2[0]);
	free(p->u2);
	free(h);
}

void tls2(void *h, float *u1, float *u2)
/*< 2*n1*n2 ==> n1*n2 >*/
{
	int n0, n1, n2, i0, i1, i2, k1, k2;
	float a, b;
	struct tag_tls2 *p = (struct tag_tls2*) h;

	n0 = p->n0;
	n1 = p->n1;
	n2 = p->n2;

#ifdef _OPENMP
#pragma omp parallel for                    \
	schedule(dynamic,10)         \
	private(i1,i2,i0,k1,a,b)
#endif
	for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	{
		for (i0=0; i0<n0; i0++)
			p->u1[i2][i1*n0+i0] = 0.0;
		for (k1=i1-p->rc1; k1<=i1+p->rc1; k1++)
		{
			if(k1<0 || k1>=n1) continue;
			a = u2[i2*n1+k1];
			b = u1[i2*n1+k1];
			p->u1[i2][i1*n0] += a*a;
			p->u1[i2][i1*n0+1] += a*b;
			if(n0 == 3)
				p->u1[i2][i1*n0+2] += b*b;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for                    \
	schedule(dynamic,n1)         \
	private(i1,i2,k2)
#endif
	for (i1=0; i1<n1*n0; i1++)
	for (i2=0; i2<n2; i2++)
	{
		p->u2[i2][i1] = 0.0;
		for (k2=i2-p->rc2; k2<=i2+p->rc2; k2++)
		{
			if(k2<0 || k2>=n2) continue;
			p->u2[i2][i1] += p->u1[k2][i1];
		}
	}

	if(p->rc3 > 0)
		recursion(p->h, p->u2[0]);

#ifdef _OPENMP
#pragma omp parallel for                    \
	schedule(dynamic,10)         \
	private(i1,i2,a,b)
#endif
	for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	{
		b = p->u2[i2][i1*n0+1];
		if(n0 == 3)
		{
			a = (p->u2[i2][i1*n0] - p->u2[i2][i1*n0+2])/2.0;
			a = sqrt(a*a+b*b) + a;
		}else{
			a = p->u2[i2][i1*n0];
		}
		if(p->ang) u1[i2*n1+i1] = atan2(b, a);
		else u1[i2*n1+i1] = b/a;
	}
}

