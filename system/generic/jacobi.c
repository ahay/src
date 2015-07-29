/* Jacobi rotation of a symmetric matrix */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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
#include <float.h>

#include <rsf.h>

static int n;
static float *aj, *ak;

void jacobi_init(int n1 /* matrix size */)
/*< initialize >*/
{
    n = n1;
    aj = sf_floatalloc(n);
    ak = sf_floatalloc(n);
}

void jacobi_close(void)
/*< free allocated storage >*/
{
    free(aj);
    free(ak);
}

float jacobi(float** a     /* matrix to rotate */, 
	     int j , int k /* indeces (k > j) */,
             float ** v    /* eigenvalues */) 
/*< Jacobi rotation >*/
{
    int i;
    float t, c, s, c2, s2, cs;

    if (k <= j) sf_error("%s: need k > j",__FILE__);
    if (fabsf(a[j][k]) < FLT_EPSILON) return 0.0f;

    t = (a[j][j]-a[k][k])/a[j][k];
    t = SF_SIG(t)/(fabsf(t)+hypotf(1.0f,t));
    c = 1.0f/hypotf(1.0f,t);
    s = t*c;
    c2=c*c;
    s2=s*s;
    cs=c*s;

    aj[j] = a[j][j]*c2+2*a[j][k]*cs+a[k][k]*s2;
    ak[k] = a[j][j]*s2-2*a[j][k]*cs+a[k][k]*c2;
    ak[j] = (a[k][k]-a[j][j])*cs+a[j][k]*(c2-s2);
    aj[k]=ak[j];

    for (i=0; i < j; i++) {
	aj[i] = a[i][k]*s+a[i][j]*c;
	ak[i] = a[i][k]*c-a[i][j]*s;
    }
    for (i=j+1; i < k; i++) {
	aj[i] = a[i][k]*s+a[j][i]*c;
	ak[i] = a[i][k]*c-a[j][i]*s;
    }
    for (i=k+1; i < n; i++) {
	aj[i] = a[k][i]*s+a[j][i]*c;
	ak[i] = a[k][i]*c-a[j][i]*s;
    }
    for (i=0; i < n; i++) {
	if (i < j) {
	    a[i][j] = aj[i];
	} else {
	    a[j][i] = aj[i];
	} 
    }
    for (i=0; i < n; i++) {
	if (i < k) {
	    a[i][k] = ak[i];
	} else {
	    a[k][i] = ak[i];
	}
    }

    if (NULL != v) {
	for (i=0; i < n; i++) {
	    aj[i] = v[j][i]*c+v[k][i]*s;
	    ak[i] = v[k][i]*c-v[j][i]*s;
	}
	for (i=0; i < n; i++) {
	    v[j][i] = aj[i];
	    v[k][i] = ak[i];
	}
    }

    return s2;
}
