/* linear phase filter bank */

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
#include "lphase.h"
#include "dsp.h"

static int nf;
static float **c;


void fbank_init(int mf, int interp)
/*< initialize >*/
{
	nf = mf;
	c = lphase(mf, interp);
}

void fbank_init2(int mf, float **a)
/*< initialize >*/
{
	nf = mf;
	c = sf_floatalloc2(2*nf+1, 2*nf+1);
	memcpy(c[0], a[0], (2*nf+1)*(2*nf+1)*sizeof(float));
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

	for(j1=0; j1<=2*nf; j1++)
	{
		firs(-nf, nf, c[j1]+nf, in, 1, n1, out[j1], 1);
	}
}


void fbank2(int n1, int n2, float**in, float ****out)
/*< 2d filter bank >*/
{
	int i1, i2, j1, j2;
	float **u1;

	u1 = sf_floatalloc2(n1, n2);

	for(j2=0; j2<=2*nf; j2++)
	for(j1=0; j1<=2*nf; j1++)
	{
		for(i1=0; i1<n1; i1++)
			firs(-nf, nf, c[j2]+nf, in[0]+i1, n1, n2, u1[0]+i1, n1);
		for(i2=0; i2<n2; i2++)
			firs(-nf, nf, c[j1]+nf, u1[i2], 1, n1, out[j2][j1][i2], 1);
	}

	free(u1[0]);
	free(u1);
}


