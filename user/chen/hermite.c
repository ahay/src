/* Hermite polynomials */

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

static int n;
static int **c;

void hermite_init(int m)
/*< initialize >*/
{
	int i1, i2;
	n = m;
	c = sf_intalloc2(n+1, n+1);

	c[0][0] = 1;
	c[1][0] = 0;	c[1][1] = 1;
	for(i2=2; i2<=m; i2++)
	{	
		for(i1=0; i1<i2; i1++)
			c[i2][i1+1] = c[i2-1][i1];
		for(i1=0; i1<i2-1; i1++)
			c[i2][i1] -= (i2-1)*c[i2-2][i1];
	}
}

void hermite_close()
/*< release the memory >*/
{
	free(c[0]);
	free(c);
}


int* hermite_coef(int i)
/*< coefficients >*/
{
	return c[i];
}

float hermite_val(int k, float x)
/*< Polynomial value >*/
{
	int i;
	float r;
	r = c[k][k];
	for(i=k-1; i>=0; i--)
		r = x*r + c[k][i];
	return r;
}



