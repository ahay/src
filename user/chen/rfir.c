/* recursive FIR filter for array signals */

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

struct tag_rfir
{
	float *c;
	float *buf, **b;
	int nf, n1;
};


void* rfir_init(int mf, float *c, int m1)
/*< FIR filter initialize >*/
{
	int i;
	struct tag_rfir *p;
	p=sf_alloc(1, sizeof(struct tag_rfir));

	p->n1 = m1;
	p->nf = mf;
	p->c = sf_floatalloc(mf);
	for(i=0; i<mf; i++) p->c[i] = c[i];

	p->b = sf_floatalloc2(m1, mf);
	p->buf = p->b[0];
	for(i=0; i<m1*mf; i++) p->buf[i] = 0.0;
	return p;
}

void rfir_close(void*h)
/*< recursive fir filter >*/
{
	struct tag_rfir *p;

	p=(struct tag_rfir*)h;
	free(p->b);
	free(p->c);
	free(p->buf);
	free(p);
}

void rfir(void *h, float*d)
/*< recursive FIR filter >*/
{
	struct tag_rfir *p;
	int i1, i2, n1, n2;
	float *pp;
	p=(struct tag_rfir*)h;

	n1 = p->n1;	
	n2 = p->nf;	

	pp = p->b[n2-1];
	for(i2=n2-1; i2>0; i2--) p->b[i2] = p->b[i2-1];
	p->b[0] = pp;
	for(i1=0; i1<n1; i1++) p->b[0][i1] = d[i1];

	for(i1=0; i1<n1; i1++)
	{
		d[i1] = 0.0;
		for(i2=0; i2<n2; i2++)
			d[i1] += p->c[i2]*p->b[i2][i1];
	}
}

