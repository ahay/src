/* filter along a dip direction */

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
#include "sinterp.h"
#include "sfilt.h"

struct tag_dipflt
{
	int nf;
	int n1;
	int n2;
	float *v1;
	sfilt filt;
	sinterp intp;
};


void *dipflt_init(int mf, int m1, int m2, 
	 char *flt, char *interp)
/*< initialize >*/
{
	struct tag_dipflt *p;
	
	p = (struct tag_dipflt*) sf_alloc(1, sizeof(struct tag_dipflt));

	p->nf = mf;
	p->n1 = m1;
	p->n2 = m2;

	p->v1 = sf_floatalloc(2*mf+1);
	p->intp = sinterp_c2f(interp);
	p->filt = sfilt_c2f(flt);
	return p;
}

void dipflt_close(void *h)
/*< release memory >*/
{
	struct tag_dipflt *p;
	p = (struct tag_dipflt*) h;
	free(p->v1);
	free(p);
}

void dipflt(void* h, float **dip, float **in, float **out)
/*< dip filter >*/
{
	int i1, i2, j1;
	struct tag_dipflt *p;
	float *pv;
	
	p = (struct tag_dipflt*) h;
	pv = p->v1+p->nf;

	for(i2=p->nf; i2<p->n2-p->nf; i2++)
	{
		for(i1=0; i1<p->n1; i1++)
		{
			for(j1=-p->nf; j1<=p->nf; j1++)
				pv[j1] = p->intp(in[i2+j1], dip[i2][i1]*j1+i1, p->n1);
			out[i2][i1] = p->filt(2*p->nf+1, p->v1);
		}
	}
}


