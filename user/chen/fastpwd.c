/* Analytical dip. */

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

static int n1, n2;
static float *a, *b, *c;

static void quadratic(float aa, float bb, float cc, 
		      float *num, float *den);

void fastpwd_init(int m1, int m2 /* data dimensions */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;

    a = sf_floatalloc(n1);
    b = sf_floatalloc(n1);
    c = sf_floatalloc(n1);
}

void fastpwd_close(void)
/*< free allocated storage >*/
{
    free(a);
    free(b);
    free(c);
}

void fastpwd(float **slice   /* [n2][n1] input data slice */,
	     float **numer   /* [n2][n1] output numerator */,
	     float **denom   /* [n2][n1] output denominator */)
/*< compute numerator and denominator for division >*/
{
    int i1, i2;
    float an, bn, cn;

    for (i1=0; i1 < n1; i1++) {
	numer[n2-1][i1] = denom[n2-1][i1] = 0.0;
	if (i1 > 0 && i1 < n1-1) {
	    a[i1] = slice[0][i1-1] - 2.0*slice[0][i1] + slice[0][i1+1];
	    b[i1] = slice[0][i1+1] - slice[0][i1-1];
	    c[i1] = slice[0][i1-1] + 4.0*slice[0][i1] + slice[0][i1+1];
	}
    }
    
    for (i2=1; i2 < n2; i2++) {
	numer[i2-1][0]    = denom[i2-1][0]    = 0.0;
	numer[i2-1][n1-1] = denom[i2-1][n1-1] = 0.0;

	for (i1=1; i1 < n1-1; i1++) {
	    an = slice[i2][i1-1] - 2.0*slice[i2][i1] + slice[i2][i1+1];
	    bn = slice[i2][i1+1] - slice[i2][i1-1];
	    cn = slice[i2][i1-1] + 4.0*slice[i2][i1] + slice[i2][i1+1];

	    quadratic((an-a[i1]),1.5*(bn+b[i1]),2.0*(cn-c[i1]),numer[i2-1]+i1,denom[i2-1]+i1);

	    a[i1] = an;
	    b[i1] = bn;
	    c[i1] = cn;
	}
    }
}

static void quadratic(float aa, float bb, float cc, 
		      float *num, float *den)
/* solution of a*x^2+2*b*x+c=0 is num/den */
{
    float d;

    d = bb*bb-aa*cc;
    
    if (d < 0.) {
	*num=0.0;
	*den=0.0;
	return;
    }

    d = sqrtf(d);

    if (bb > 0) {
	*den = bb+d;
    } else {
	*den = bb-d;
    }
	
	*num = -cc;
}
