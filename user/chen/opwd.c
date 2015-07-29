/* circle-interpolating pwd*/

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

static int nf1, nf2;
static float *b1, *b2;
static float r, **c;

void opwd_init(int mf1, int mf2, char *interp, float rad)
/*< initialize >*/
{

	nf1 = mf1;
	nf2 = mf2;
	r = rad;
	b1 = sf_floatalloc(mf1+mf2+1);
	b2 = sf_floatalloc(mf1+mf2+1);

	c = lphpoly(mf1, mf2, interp);
}

void opwd(int n1, int n2, float **in, sf_complex **p, float **out)
/*< apply circle interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=nf1; i2<n2-nf2; i2++)
	for(i1=nf1; i1<n1-nf2; i1++)
	{
		lphpoly_coef(nf1+nf2, c, r*cimagf(p[i2][i1]), b1, nf1+nf2, false);
		lphpoly_coef(nf1+nf2, c, r*crealf(p[i2][i1]), b2, nf1+nf2, false);
		c1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, 1, n1);
		c2 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, -1, -n1);
		out[i2][i1] = c2 - c1;
	}
}

void opwdpd(int n1, int n2, float **in, 
	sf_complex **p, float **out, int id)
/*< partial derivative filter of circle interpolating PWD >*/
{
	int i1, i2;
	float p1, p2, c1, z1;

	for(i2=nf1; i2<n2-nf2; i2++)
	for(i1=nf1; i1<n1-nf2; i1++)
	{
		p1 = r*cimagf(p[i2][i1]);
		p2 = r*crealf(p[i2][i1]);
		if(id==0)
		{
			lphpoly_coef(nf1+nf2, c, p1, b1, nf1+nf2, true);
			lphpoly_coef(nf1+nf2, c, p2, b2, nf1+nf2, false);
		} else{
			lphpoly_coef(nf1+nf2, c, p1, b1, nf1+nf2, false);
			lphpoly_coef(nf1+nf2, c, p2, b2, nf1+nf2, true);
		}
		z1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, 1, n1);
		c1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, -1, -n1);
		out[i2][i1] = (c1 - z1);/*(id==0?p2:-p1);*/
	}
}



void opwd_freq(sf_complex dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2;
	sf_complex c1, c2, z1, z2;

	lphpoly_coef(nf1+nf2, c, r*cimagf(dip), b1, nf1+nf2, false);
	lphpoly_coef(nf1+nf2, c, r*crealf(dip), b2, nf1+nf2, false);

	for(i2=-nk; i2<=nk; i2++)
	{
		c2 = fir_freq(-nf1, nf2, b2+nf1, 0.5*i2/nk);
		for(i1=-nk; i1<=nk; i1++)
		{
			c1 = fir_freq(-nf1, nf2, b1+nf1, 0.5*i1/nk);
			z1 = conjf(c1);	z2 = conjf(c2);
#ifdef SF_HAS_COMPLEX_H
			if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/(z1*z2);
			else	out[i2+nk][i1+nk] = z1*z2 - c1*c2;
#else
			if(iir)	out[i2+nk][i1+nk] = sf_cadd(sf_cmplx(1.0,0.0),sf_cneg(sf_cdiv(sf_cmul(c1,c2),sf_cmul(z1,z2))));
			else	out[i2+nk][i1+nk] = sf_cadd(sf_cmul(z1,z2),sf_cneg(sf_cmul(c1,c2)));
#endif
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

