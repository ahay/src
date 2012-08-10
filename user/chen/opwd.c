/* circle-interpolating pwd*/

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
static float *b1, *b2;
static float r, **c;

void opwd_init(int interp, int mf, float rad)
/*< initialize >*/
{

	itp = interp;
	nf = mf;
	r = rad;
	b1 = sf_floatalloc(2*mf+1);
	b2 = sf_floatalloc(2*mf+1);

	c = lphase(mf, interp);
}

void opwd(int n1, int n2, float **in, sf_complex **p, float **out)
/*< apply circle interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		lphase_filt(nf, c, r*cimagf(p[i2][i1]), b1, 2*nf, false);
		lphase_filt(nf, c, r*crealf(p[i2][i1]), b2, 2*nf, false);
		c1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, 1, n1);
		switch(itp)
		{
		case 1:
			c2 = in[i2][i1];
			break;
		case 2:
			c2 = fir2(-nf, nf, c[0]+nf, -nf, nf, c[0]+nf, in[i2]+i1, 1, n1);
			break;
		default:
			c2 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, -1, -n1);
		}
		out[i2][i1] = c2 - c1;
	}
}

void opwdpd(int n1, int n2, float **in, 
	sf_complex **p, float **out, int id)
/*< partial derivative filter of circle interpolating PWD >*/
{
	int i1, i2;
	float p1, p2, c1, z1;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		p1 = r*cimagf(p[i2][i1]);
		p2 = r*crealf(p[i2][i1]);
		if(id==0)
		{
			lphase_filt(nf, c, p1, b1, 2*nf, true);
			lphase_filt(nf, c, p2, b2, 2*nf, false);
		} else{
			lphase_filt(nf, c, p1, b1, 2*nf, false);
			lphase_filt(nf, c, p2, b2, 2*nf, true);
		}
		z1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, 1, n1);
		switch(itp)
		{
		case 1:
			c1 = 0.0;
			break;
		case 2:
			c1 = 0.0;
			break;
		default:
			c1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, -1, -n1);
		}
		out[i2][i1] = (c1 - z1);//*(id==0?p2:-p1);
	}
}



void opwd_freq(sf_complex dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2;
	sf_complex c1, c2, z1, z2;

	lphase_filt(nf, c, r*cimagf(dip), b1, 2*nf, false);
	lphase_filt(nf, c, r*crealf(dip), b2, 2*nf, false);

	for(i2=-nk; i2<=nk; i2++)
	{
		c2 = fir_freq(-nf, nf, b2+nf, 0.5*i2/nk);
		for(i1=-nk; i1<=nk; i1++)
		{
			c1 = fir_freq(-nf, nf, b1+nf, 0.5*i1/nk);
			switch(itp)
			{
			case 1:
				z1 = 1.0; z2 = 1.0;
				break;
			case 2:
				z1 = fir_freq(-nf, nf, c[0]+nf, 0.5*i1/nk);
				z2 = fir_freq(-nf, nf, c[0]+nf, 0.5*i2/nk);
				break;
			default:
				z1 = conj(c1);	z2 = conj(c2);
			}
			if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/(z1*z2);
			else	out[i2+nk][i1+nk] = z1*z2 - c1*c2;
		}
	}

}



void opwd_close()
/*< release memory >*/
{
	free(b1);
	free(b2);
	free(c[0]);
	free(c);
}

