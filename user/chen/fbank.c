/* linear phase filter bank */

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

static int m, n;
static float **c;


void fbank_init(int mf, int nf, char* interp)
/*< initialize >*/
{
	m = mf;
	n = nf;
	c = lphpoly(m, n, interp);
}

void fbank_init2(int mf, int nf, float **a)
/*< initialize >*/
{
	m = mf;
	n = nf;
	c = sf_floatalloc2(m+n+1, m+n+1);
	memcpy(c[0], a[0], (m+n+1)*(m+n+1)*sizeof(float));
}

void fbank_close()
/*< release memory >*/
{
	free(c[0]);
	free(c);
}


void fbank1(int n1, float *in, float **out)
/*< 1d filter bank >*/
{
	int j1;

	for(j1=0; j1<=(m+n); j1++)
	{
		firs(-m, n, c[j1]+m, in, 1, n1, out[j1], 1);
	}
}


void fbank2(int n1, int n2, float**in, float ****out)
/*< 2d filter bank >*/
{
	int i1, i2, j1, j2;
	float **u1;

	u1 = sf_floatalloc2(n1, n2);

	for(j2=0; j2<=(m+n); j2++)
	for(j1=0; j1<=(m+n); j1++)
	{
		for(i1=0; i1<n1; i1++)
			firs(-m, n, c[j2]+m, in[0]+i1, n1, n2, u1[0]+i1, n1);
		for(i2=0; i2<n2; i2++)
			firs(-m, n, c[j1]+m, u1[i2], 1, n1, out[j2][j1][i2], 1);
	}

	free(u1[0]);
	free(u1);
}


