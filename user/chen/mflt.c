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

struct tag_mflt
{
	float *bi;
	int n1, nw;
	int *pi;
};


void* mflt_init(int m1, int mw)
/*< median filter initialize >*/
{
	struct tag_mflt *p;
	p=sf_alloc(1, sizeof(struct tag_mflt));

	p->n1 = m1;
	p->nw = mw;
	p->bi = sf_floatalloc(m1+2*mw);
	p->pi = sf_intalloc(2*mw+1);
	return p;
}

void mflt_close(void*h)
/*< recursive median filter >*/
{
	struct tag_mflt *p;

	p=(struct tag_mflt*)h;
	free(p->pi);
	free(p->bi);
	free(p);
}


void mflt(void *h, float*d, int d1)
/*< recursive median filter >*/
{
	struct tag_mflt *p;
	int i1, i2, j1, temp, chg, *pi, n1, nw;
	float *pb;
	
	p=(struct tag_mflt*)h;
	nw=p->nw;
	n1=p->n1;
	pi=p->pi+nw;
	pb=p->bi+nw;

	for(i1=0; i1<n1; i1++)
		pb[i1] = d[i1*d1];
	for(i1=1; i1<=nw; i1++)
	{
		pb[0-i1] = pb[0];
		pb[n1-1+i1] = pb[n1-1];
		pi[-i1] = -i1;
		pi[i1] = i1;
	}
	pi[0] = 0;

	/* initialiaze */
	for(j1=nw; j1>=-nw; j1--)
	{
		chg=0;
		for(i1=-nw; i1<j1; i1++)
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

	/* update */
	for(i2=1; i2<n1; i2++)
	{
		pb++;
		for(i1=-nw; i1<=nw; i1++)
		{
			pi[i1] -= 1;
			if(pi[i1] < -nw)
			{
				j1 = i1;
				pi[i1] = nw;
			}
		}	
		for(i1=j1; i1<nw; i1++)
		{
			if(pb[pi[i1]] > pb[pi[i1+1]])
			{
				temp = pi[i1];
				pi[i1] = pi[i1+1];
				pi[i1+1] = temp;
			}else break;
		}
		for(i1=j1; i1>-nw; i1--)
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

