/* general recursive operator */

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

typedef void (*op_func)(float *out, float **in, int n1, int n2, int*par);
/*^*/

struct tag_recursion
{
	float **b, *buf;
	int n1, n2;
	op_func op;
	int *par;
};

void* recursion_init(int n1, int n2, op_func pop, int *para)
/*< initialize >*/
{
	int i;
	struct tag_recursion *p;
	p=sf_alloc(1, sizeof(struct tag_recursion));

	p->n1 = n1; /* array size */
	p->n2 = n2; /* operator size */
	p->op = pop;

	p->b = sf_floatalloc2(n1, p->n2);
	p->buf = p->b[0];
	p->par = para;
	for(i=0; i<n1*p->n2; i++) p->buf[i] = 0.0;
	return p;
}

void recursion_push(void *h, float *d)
/*< push data buf >*/
{
	struct tag_recursion *p;
	int i1, i2;
	float *pp;
	p=(struct tag_recursion*)h;
	
	/* update */
	pp = p->b[p->n2-1];
	for(i2=p->n2-1; i2>0; i2--) p->b[i2] = p->b[i2-1];
	p->b[0] = pp;
	for(i1=0; i1<p->n1; i1++) p->b[0][i1] = d[i1];
}

void recursion_close(void*h)
/*< release recursive operator >*/
{
	struct tag_recursion *p;

	p=(struct tag_recursion*)h;
	free(p->b);
	free(p->buf);
	free(p);
}

void recursion(void *h, float*d)
/*< recursive operator >*/
{
	struct tag_recursion *p;
	p=(struct tag_recursion*)h;
	recursion_push(h, d);
/*	p->upd(p->b, p->n2); */
	p->op(d, p->b, p->n1, p->n2, p->par);
}

