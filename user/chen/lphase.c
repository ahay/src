/* Linear phase filters */

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
#include "poly.h"

float **lphlag(int n)
/*< Linear phase filter by noncausal Lagrange approximation
[    1     p    p^2  |
   ------------------+---
    0.0   0.5   0.5  |  Z
    1.0   0.0  -1.0	 |
    0.0  -0.5   0.5  |  Z^{-1} ]
 >*/
{
	float **p, a;
	int i1, i2, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	for(i1=-n; i1<=n; i1++)
	{
		p[n-i1][m] = 1.0;
		for(i2=-n; i2<i1; i2++)
		{
			p[n-i1][m] /= (i1-i2);
			p[n-i1][i2+n] = i2;
		}
		for(i2=i1+1; i2<=n; i2++)
		{
			p[n-i1][m] /= (i1-i2);
			p[n-i1][i2+n-1] = i2;
		}
		poly(m, p[n-i1]);
	}
	for(i2=0; i2<=m; i2++)
	for(i1=0; i1<i2; i1++)
	{
		a = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = a;
	}
	return p;
}


float **lphbsp(int n)
/*< Linear phase filter by noncausal B-Spline approximation
[    1     p    p^2  |
   ------------------+---
   0.25   1.0   1.0  |  Z
    1.5   0.0  -2.0	 |
   0.25  -1.0   1.0  |  Z^{-1} ]
 >*/
{
	float **p, *a, *b, c;
	int i1, i2, j, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	a = sf_floatalloc(m+1);
	b = sf_floatalloc(m+1);
	
	a[0] = 1.0;
	b[0] = 1.0;
	for(i2=1; i2<=m; i2++)
	{
		a[i2] =-a[i2-1] * (m+2-i2) / i2;
		b[i2] = b[i2-1] * (m+1-i2) / i2;
	}

	for(i1=-n; i1<=n; i1++)
	{
		for(i2=0; i2<=m; i2++)
		{
			p[n-i1][i2] = 0.0;
			for(j=0; j<=n-i1; j++)
				p[n-i1][i2] += a[j]*powf(n+0.5-i1-j, m-i2);
			p[n-i1][i2] *= b[i2];
		}
	}
	for(i2=0; i2<=m; i2++)
	for(i1=0; i1<i2; i1++)
	{
		c = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = c;
	}
	free(a);
	free(b);
	return p;
}


float **lphmf(int n)
/*< Linear phase filter by noncausal maxflat approximation 
           B(Z)
  H(Z) = --------
          B(1/Z)
[    1     p    p^2  |
   ------------------+---
    1.0   1.5   0.5  |  Z
    4.0   0.0  -1.0	 |
    1.0  -1.5   0.5  |  Z^{-1} ]
 >*/
{
	float **p, a;
	int i1, i2, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	for(i1=-n; i1<=n; i1++)
	{
		p[i1+n][m] = ((n+i1)%2==0?1.0:-1.0);
		for(i2=0; i2<i1+n; i2++)
		{
			p[i1+n][m] /= (i2+1);
			p[i1+n][i2] = -(i2-m);
		}
		for(i2=i1+n+1; i2<=m; i2++)
		{
			p[i1+n][m] /= (i2-i1-n);
			p[i1+n][i2-1] = -i2;
		}
		poly(m, p[i1+n]);
	}
	for(i2=0; i2<=m; i2++)
	for(i1=0; i1<i2; i1++)
	{
		a = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = a;
	}

	for(i1=1, a=1.0; i1<=m; i1++)
		a = a*i1/(m+i1);

	for(i2=0; i2<=m; i2++)
	for(i1=0; i1<=m; i1++)
		p[i1][i2] *= a;

	return p;
}

float **lphase(int n, int interp)
/*< linear phase filter bank >*/
{
	switch(interp)
	{
	case 1:
		return lphlag(n);
	case 2:
		return lphbsp(n);
	default:
		return lphmf(n);
	}
}


void lphase_filt(int n, float **c, float delay, 
	float *out, int np, bool der)
/*< filter >*/
{
	int i1, i2;
	double b;

	for(i1=0; i1<=2*n; i1++)
	{
		b = der?0:c[0][i1];
		for(i2=1; i2<=np; i2++)
			b += c[i2][i1] * 
				(der?(powf(delay, i2-1)*i2):powf(delay, i2));
		out[i1] = b;
	}
}


