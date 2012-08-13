/* filter along a dip direction */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

struct tag_dipfilt
{
	int nf;
	int n1;
	int n2;
	float *v1;
	sfilt filt;
	sinterp interp1;
	sinterp interp2;
};


void *dipfilt_init(int mf, int m1, int m2, 
	 sfilt flt, sinterp int1, sinterp int2)
/*< initialize >*/
{
	struct tag_dipfilt *p;
	
	p = (struct tag_dipfilt*) sf_alloc(sizeof(struct tag_dipfilt), 1);

	p->nf = mf;
	p->n1 = m1;
	p->n2 = m2;

	p->v1 = sf_floatalloc(2*mf+1);
	p->filt = flt;
	p->interp1 = int1;
	p->interp2 = int2;
	return p;
}


void dipfilt(void* h, float **dip, float **in, float **out)
/*< initialize >*/
{
	int i1, i2, i3;
	struct tag_dipfilt *p;
	float *pv, d1;
	
	pv = p->v1+p->nf;
	p = (struct tag_dipfilt*) h;

	for(i2=p->nf; i2<p->n2-p->nf; i2++)
	{
		for(i1=p->nf; i1<p->n1-p->nf; i1++)
		{
			pv[0] = in[i2][i1];
			d1 = dip[i2][i1];
			for(i3=1; i3<=p->nf; i3++)
			{
				pv[i3] = p->interp1(in[i2+i3], i1+d1, p->n1);
				d1 = p->interp2(dip[i2+i3], i1+d1, p->n1);
			}
			d1 = -dip[i2][i1];
			for(i3=1; i3<=p->nf; i3++)
			{
				pv[-i3] = p->interp1(in[i2-i3], i1+d1, p->n1);
				d1 = -p->interp2(dip[i2+i3], i1+d1, p->n1);
			}
			out[i2][i1] = p->filt(2*p->nf+1, p->v1);
		}
	}
}


