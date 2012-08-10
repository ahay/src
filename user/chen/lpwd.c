/* line-interpolating pwd*/

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

static int nf, itp;
static int n1, n2;
static float **c, *b1;

void lpwd_init(int interp, int mf, int m1, int m2)
/*< initialize >*/
{

	itp = interp;
	nf = mf;
	n1 = m1;
	n2 = m2;

	c = lphase(nf, interp);

	b1 = sf_floatalloc(2*mf+1);
}

void lpwc(float **in, float **out, float **p)
/*< apply line interpolating PWC >*/
{
/*
	int i1, i2, j1;

	for(i2=0; i2<n2; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		lphase_filt(nf, c, p[i2][i1], b1, 2*nf, false);
		out[i2][i1] = 0.0;
		for(j1=-nf; j1<=nf; j1++)
			out[i2][i1] += b1[j1+nf]*in[i2][i1-j1];
	}
*/
}

void lpwd(float **in, float **out, float **p)
/*< apply line interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=0; i2<n2-1; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		lphase_filt(nf, c, p[i2][i1], b1, 2*nf, false);
		switch(itp)
		{
		case 1:
			c1  = in[i2][i1];
			break;
		case 2:
			c1  = fir(-nf, nf, c[0]+nf, in[i2]+i1, 1);
			break;
		default:
			c1  = fir(-nf, nf, b1+nf, in[i2]+i1, -1);
		}
		c2 = fir(-nf, nf, b1+nf, in[i2+1]+i1, 1);
		out[i2][i1]  = c1 - c2;
	}
}


void lpwdd(float **in, float **out, float **p)
/*< derivative filter of line interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=0; i2<n2-1; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		lphase_filt(nf, c, p[i2][i1], b1, 2*nf, true);
		switch(itp)
		{
		case 1:
		case 2:
			c1  = 0.0;
			break;
		default:
			c1  = fir(-nf, nf, b1+nf, in[i2]+i1, -1);
		}
		c2 = fir(-nf, nf, b1+nf, in[i2+1]+i1, 1);
		out[i2][i1]  = c1 - c2;
	}
}

void lpwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of line-interpolating PWD >*/
{
	int i1, i2;
	sf_complex c1, c2, z1;

	lphase_filt(nf, c, tan(dip), b1, 2*nf, false);

	for(i2=-nk; i2<=nk; i2++)
	{
		c2 = cexpf(sf_cmplx(0., SF_PI*i2/nk));
		for(i1=-nk; i1<=nk; i1++)
		{
			c1 = fir_freq(-nf, nf, b1+nf, 0.5*i1/nk);
			switch(itp)
			{
			case 1:
				out[i2+nk][i1+nk] = 1.0 - c2*c1;
				break;
			case 2:
				z1 = fir_freq(-nf, nf, c[0]+nf, 0.5*i1/nk);
				if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/z1;
				else	out[i2+nk][i1+nk] = z1 - c1*c2;
				break;
			default:
				z1 = conj(c1);
				if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/z1;
				else	out[i2+nk][i1+nk] = z1 - c1*c2;
			}
		}
	}
}


void lpwd_close()
/*< release memory >*/
{
	free(b1);
	free(c[0]);
	free(c);
}

