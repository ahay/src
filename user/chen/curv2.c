/* Curv2 */

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
#include "lphpoly.h"
#include "dsp.h"

static int m, n, n1, n2, iter;
static float **ut, **ux, **u2x2, **u2tx;
static float *b1, *b2, **c;

void curv2_init(int mf, int nf, char* interp,
	int m1, int m2, int niter)
/*< initialize >*/
{
	int i1;

	m = mf;
	n = nf;
	c = lphpoly(m, n, interp);
	b1 = c[1]+m;
	b2 = c[2]+m;
	for(i1=-m; i1<=n; i1++)
		b2[i1] /= 2;
	n1 = m1;
	n2 = m2;
	iter = niter;
	ut  = sf_floatalloc2(n1, n2);
	ux  = sf_floatalloc2(n1, n2);
	u2x2  = sf_floatalloc2(n1, n2);
	u2tx  = sf_floatalloc2(n1, n2);

}

void curv2_close()
/*< release memory >*/
{
	free(*ut);
	free(ut);
	free(*ux);
	free(ux);
	free(*u2tx);
	free(u2tx);
	free(*u2x2);
	free(u2x2);
	free(*c);
	free(c);
}

void curv2(float **u, float **zx)
/*< curv2 >*/
{
	int i1, i2, i3;
	for(i1=0; i1<n1; i1++)
	{
		firs(m,n,b1,u[0]+i1,n1,n2,ux[0]+i1,n1);
		firs(m,n,b2,u[0]+i1,n1,n2,u2x2[0]+i1,n1);
	}
	for(i2=0; i2<n2; i2++)
	{
		firs(m,n,b1,u[i2],1,n1,ut[i2],1);
		firs(m,n,b1,ux[i2],1,n1,u2tx[i2],1);
		for(i1=0; i1<n1; i1++) u[i2][i1] = 0.0;
	}

	for(i3=0; i3<iter; i3++)	// iterations
	{
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
		{
			zx[i2][i1] = -(ut[i2][i1]*ux[i2][i1] + u2tx[i2][i1]*u2x2[i2][i1]);
			zx[i2][i1] /= (ut[i2][i1]*ut[i2][i1] + u2tx[i2][i1]*u2tx[i2][i1]);
		}
		for(i1=0; i1<n1; i1++)
		firs(m,n,b1,zx[0]+i1,n1,n2,u[0]+i1,n1);
	}
}



