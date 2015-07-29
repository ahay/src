/* Linear PHase approximation by POLYnomial coefficients filters */

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
#include "poly.h"

float **lphlag(int m, int n)
/*< Linear phase filter by noncausal Lagrange approximation
[    1     p    p^2  |
   ------------------+---
    0.0   0.5   0.5  |  Z^{m}
    1.0   0.0  -1.0	 |
    0.0  -0.5   0.5  |  Z^{-n} ]
 >*/
{
	float **p, c;
	int i, k, n1, i1, i2;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	for(k=-m; k<=n; k++)
	{
		p[k+m][n1] = 1.0;
		for(i=-m; i<k; i++)
		{
			p[k+m][n1] /= (i-k);
			p[k+m][i+m] = -i;
		}
		for(i=k+1; i<=n; i++)
		{
			p[k+m][n1] /= (i-k);
			p[k+m][i+m-1] = -i;
		}
		poly(n1, p[k+m]);
	}

	for(i2=0; i2<=n1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		c = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = c;
	}

	return p;
}

static double factorial(int s, int e)
{
	double c=1.0;
	int i;
	for(i=s; i<=e; i++) c*=i;
	return c;	
}

float **lphbsp(int m, int n)
/*< Linear phase filter by noncausal B-Spline approximation
[    1     p    p^2  |
   ------------------+---
   0.125   1.0   1.0  |  Z
   0.75    0.0  -2.0	 |
   0.125  -1.0   1.0  |  Z^{-1} ]
 >*/
{
	float **p;
	double c, c0, c1, c2, c3;
	int i, j, k, n1;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	c0 = factorial(1, n1+1);
	for(j=0; j<=n1; j++)
	{
		c1 = factorial(1, j) * factorial(1, n1-j);
		for(k=-m; k<=n; k++)
		{
			c = 0;
			for(i=0; i<k+1+0.5*n1; i++)
			{
				c2 = (i%2 == 1)? -1.0 : 1.0;
				c2 *= factorial(1, i) * factorial(1, n1+1-i);
				c3 = pow((n1+1.0)/2+k-i, n1-j);
				c += c0*c3/(c1*c2);
			}
			p[j][k+m] = c;	
		}
	}
	return p;
}


float **lphmf(int m, int n)
/*< Linear phase filter by noncausal maxflat approximation 
[    1     p    p^2  |
   ------------------+---
    0.1667   0.5   0.3333  |  Z
    0.6667   0.0  -0.6667  |
    0.1667  -0.5   0.3333  |  Z^{-1} ]
 >*/
{
	float **p;
	int i, k, n1, i1, i2;
	double c1, c2;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	c1 = factorial(1, n1);
	c1 = c1/factorial(n1+1, 2*n1);
	for(k=-m; k<=n; k++)
	{
		c2 = ((m+k)%2 == 1)? -1.0 : 1.0;
		c2 *= factorial(1, m+k) * factorial(1, n-k);
		p[k+m][n1] =c1/c2;
		for(i=0; i<m+k; i++)
			p[k+m][i] = (2*m-i);
		for(i=0; i<n-k; i++)
			p[k+m][i+m+k] = (i-2*n);
		poly(n1, p[k+m]);
	}
	for(i2=0; i2<=n1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		c1 = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = c1;
	}

	for(i1=n1; i1>=1; i1--)
	{
		for(i=n1; i>=i1; i--)	
		for(k=-m; k<=n; k++)
		p[i][k+m] *= 2.0;
	}

	return p;
}

float **lphpoly(int m, int n, char* interp)
/*< linear phase filter bank >*/
{
	if(strcmp(interp,"lagrange")==0)	return lphlag(m, n);
	else if(strcmp(interp,"bspline")==0)	return lphbsp(m, n);
	else if(strcmp(interp,"maxflat")==0)	return lphmf(m, n);
	else return NULL;
}


void lphpoly_coef(int nf, float **c, float delay, 
	float *out, int np, bool der)
/*< filter >*/
{
	int k, i;
	double b;

	for(k=0; k<=nf; k++)
	{
		b = der?0:c[0][k];
		for(i=1; i<=np; i++)
			b += c[i][k] * 
				(der?(powf(delay, i-1)*i):powf(delay, i));
		out[k] = b;
	}
}


