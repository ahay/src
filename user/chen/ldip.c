/* omnidirectional dip estimation by PWD */

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
#include "lpwd.h"

static int n1, n2, nf1, nf2;
static float **u1, **u2, **u3;
static bool verb, use_divn;

void ldip_init(char *interp, int mf1, int mf2,
	int m1, int m2, int *rect, int niter, bool vb)
/*< initialize >*/
{
	int n, nn[4];
	nf1 = mf1;
	nf2 = mf2;
	n1 = m1;
	n2 = m2;
	verb=vb;

	u1 = sf_floatalloc2(n1,n2);
	u2 = sf_floatalloc2(n1,n2);
	u3 = sf_floatalloc2(n1,n2);

	lpwd_init(nf1, nf2, n1, n2, interp);
	if(rect[0]>0 && rect[1]>0)
	{
		n = n1*n2;
		nn[0] = n1;
		nn[1] = n2;
		sf_divn_init (2, n, nn, rect, niter, vb);
		use_divn=true;
	}else use_divn=false;
}

void ldip_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(u3[0]);
	free(u3);
	lpwd_close();
	if( use_divn)
		sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+0.0001))

void ldip(float **in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double norm;

	for (it=0; it<nit; it++)
	{
		lpwd(in, u1, dip);
		lpwdd(in, u2, dip);

		if(verb)
		{
			for(i1=0, norm=0.0; i1<n1*n2; i1++)
				norm += (u1[0][i1]*u1[0][i1]);
			sf_warning("%d %g", it+1, sqrtf(norm/n1/n2));
		}
		if(use_divn)
		{
			for(i1=0, norm=0.0; i1<n1*n2; i1++)
				norm += (u2[0][i1]*u2[0][i1]);
			norm=sqrtf(norm/(n1*n2));
			for(i1=0; i1<n1*n2; i1++)
			{
				u1[0][i1] /= norm;
				u2[0][i1] /= norm;
			}

			sf_divn(u1[0], u2[0], u3[0]);
		}else{
			for(i1=0; i1<n1*n2; i1++)
				u3[0][i1] -= divn(u1[0][i1], u2[0][i1]);
		}
		for(i1=0; i1<n1*n2; i1++)
			dip[0][i1] -= eta*u3[0][i1];

	}
}



