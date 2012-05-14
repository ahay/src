/* omnidirectional plane-wave destruction */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

#include "opwd.h"

#ifndef _opwd_h

typedef void (*interpolate)(float,float*);
/*^*/

#endif

static int n1, n2;
static float **res;

void lagrange(float p, float *f)
/*< Lagrange interpolator >*/
{
    f[0]=0.5f*p*(p-1.0f);
    f[1]=1.0f-p*p;
    f[2]=0.5f*p*(p+1.0f);
}

void lagrange_der(float p, float *f)
/*< Lagrange interpolator derivative >*/
{
    f[0]=p-0.5f;
    f[1]=-2.0f*p;
    f[2]=p+0.5f;
}

void bspline(float p, float *f)
/*< B-spline interpolator >*/
{
    f[0]=0.125*(2*p-1)*(2*p-1);
    f[1]=0.75-p*p;
    f[2]=0.125*(2*p+1)*(2*p+1);
}

void opwd_init(int m1, int m2 /* data dimensions */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    res = sf_floatalloc2(n2,n1);
}

void opwd_close(void)
/*< free allocated storage >*/
{
    free(*res);
    free(res);
}

void opwd_filter(interpolate int1, interpolate int2 /* interpolators */,
		 sf_tris t1, sf_tris t2             /* prefilters */,
		 float **p1, float **p2             /* slopes [n2][n1] */,
		 float **inp, float **out           /* [n2][n1] */)
/*< filter input >*/
{
    int i1, i2;
    float *trace, filt[3];

    for (i2=0; i2 < n2; i2++) {
	res[0][i2] = res[n1-1][i2] = 0.;
	trace = inp[i2];
	if (NULL != t1) sf_tridiagonal_solve(t1,trace);
	for (i1=1; i1 < n1-1; i1++) {
	    int1(p1[i2][i1],filt);
	    res[i1][i2] = filt[0]*trace[i1-1]+filt[1]*trace[i1]+filt[2]*trace[i1+1];
	}
	out[i2][0] = out[i2][n1-1] = 0.;
    }
    for (i1=0; i1 < n1; i1++) {
	out[0][i1] = out[n2-1][i1] = 0.;
	trace = res[i1];
	if (NULL != t2) sf_tridiagonal_solve(t2,trace);
	for (i2=1; i2 < n2-1; i2++) {
	    int2(p2[i2][i1],filt);
	    out[i2][i1] = filt[0]*trace[i2-1]+filt[1]*trace[i2]+filt[2]*trace[i2+1];
	}
    }
}
