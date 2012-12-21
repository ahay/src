/* recursive median filter */

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

struct tag_rmedian
{
	float *bi;
	int n1, nw1, nw2;
	int *pi;
};


void* rmedian_init(int m1, int mw1, int mw2)
/*< median filter initialize >*/
{
	struct tag_rmedian *p;
	p=sf_alloc(1, sizeof(struct tag_rmedian));

	p->n1 = m1;
	p->nw1 = mw1;
	p->nw2 = mw2;
	p->bi = sf_floatalloc(m1+mw1+mw2);
	p->pi = sf_intalloc(mw1+mw2+1);
	return p;
}

void rmedian_close(void*h)
/*< recursive median filter >*/
{
	struct tag_rmedian *p;

	p=(struct tag_rmedian*)h;
	free(p->pi);
	free(p->bi);
	free(p);
}


void rmedian(void *h, float*d, int d1)
/*< recursive median filter >*/
{
	struct tag_rmedian *p;
	int i1, i2, j1, temp, chg, *pi, n1, nw1, nw2;
	float *pb;
	
	p=(struct tag_rmedian*)h;
	nw1=p->nw1;
	nw2=p->nw2;
	n1=p->n1;
	pi=p->pi+nw1;
	pb=p->bi+nw1;

	for(i1=0; i1<n1; i1++)
		pb[i1] = d[i1*d1];
	for(i1=1; i1<=nw1; i1++)
	{
		pb[-i1] = pb[0];
		pi[-i1] = -i1;
	}
	for(i1=1; i1<=nw2; i1++)
	{
		pb[n1-1+i1] = pb[n1-1];
		pi[i1] = i1;
	}
	pi[0] = 0;

	// initialiaze
	for(j1=nw2; j1>=-nw1; j1--)
	{
		chg=0;
		for(i1=-nw1; i1<j1; i1++)
			if(pb[pi[i1]] > pb[pi[i1+1]])
			{
				temp = pi[i1];
				pi[i1] = pi[i1+1];
				pi[i1+1] = temp;
				chg=1;
			}
		if(chg==0) break;
	}
	d[0]=pb[pi[0]];

	// update
	for(i2=1; i2<n1; i2++)
	{
		pb++;
		for(i1=-nw1; i1<=nw2; i1++)
		{
			pi[i1] -= 1;
			if(pi[i1] < -nw1)
			{
				j1 = i1;
				pi[i1] = nw2;
			}
		}	
		for(i1=j1; i1<nw2; i1++)
		{
			if(pb[pi[i1]] > pb[pi[i1+1]])
			{
				temp = pi[i1];
				pi[i1] = pi[i1+1];
				pi[i1+1] = temp;
			}else break;
		}
		for(i1=j1; i1>-nw1; i1--)
		{
			if(pb[pi[i1]] < pb[pi[i1-1]])
			{
				temp = pi[i1];
				pi[i1] = pi[i1-1];
				pi[i1-1] = temp;
			}else break;
		}
		d[i2*d1]=pb[pi[0]];
	}
}


