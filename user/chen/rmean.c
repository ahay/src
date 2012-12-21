/* recursive mean filter */

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

struct tag_rmean
{
	float *bi;
	int n1, nw1, nw2;
};


void* rmean_init(int m1, int mw1, int mw2)
/*< median filter initialize >*/
{
	struct tag_rmean *p;
	p=sf_alloc(1, sizeof(struct tag_rmean));

	p->n1 = m1;
	p->nw1 = mw1;
	p->nw2 = mw2;
	p->bi = sf_floatalloc(m1+mw1+mw2);
	return p;
}

void rmean_close(void*h)
/*< recursive median filter >*/
{
	struct tag_rmean *p;

	p=(struct tag_rmean*)h;
	free(p->bi);
	free(p);
}

void rmean(void *h, float*d, int d1)
/*< recursive mean filter >*/
{
	struct tag_rmean *p;
	int i1, i2, n1, nw1, nw2;
	float *pb, t;
	
	p=(struct tag_rmean*)h;
	nw1=p->nw1;
	nw2=p->nw2;
	n1=p->n1;
	pb=p->bi+nw1;

	for(i1=0; i1<n1; i1++) pb[i1] = d[i1*d1];
	for(i1=1; i1<=nw1; i1++) pb[-i1] = pb[0];
	for(i1=1; i1<=nw2; i1++) pb[n1-1+i1] = pb[n1-1];

	for(i2=-nw1, t=0.0; i2<=nw2; i2++)	t += pb[i2];
	d[0] = t/(nw1+nw2+1);
	for(i1=1; i1<n1; i1++)
	{
		t += pb[i1+nw2]-pb[i1-nw1];
		d[i1*d1] = t/(nw1+nw2+1);
	}
}

