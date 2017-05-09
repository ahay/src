/* cubic spline interpolation */
/*
  Copyright (C) 2011 KAUST

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

/* Reference: [Numerical mathematics and computing (6ed) by Cheney and Kincaid], equation 10 page 390 */

#undef DEBUG

#include <rsf.h>

#include "spline.h"
#ifndef _spline_h

typedef struct{
	int n;
	const float *x,*f; /* [n]   data f(x)    */
	float *a,*b,*c,*d; /* [n-1] coefficients */
	bool linear;       /* linear or cubic    */
} Spl;
/*^*/

#endif

static void mysgtsv(int n,int nrhs,float *dl,float *d,float *du,float *b,int ldb)
{
	int info;
	sgtsv_(&n,&nrhs,dl,d,du,b,&ldb,&info);
	if (info)
		sf_error("sgtsv failed with info=%d",info);
}

Spl *spline_init(const float *ff, /* [n] */
				 const float *xx, /* [n] sorted */
				 int nn,
				 bool linear)
/*< initialize >*/
{
	int i,n;
	float tmp,*h,*b,*u,*v,*A,*B,*C,*D,*f,*x;
	Spl *spl;

	spl = sf_alloc(1,sizeof(*spl));

	/* avoid dx=0 */
	f = sf_floatalloc(nn); f[0] = ff[0];
	x = sf_floatalloc(nn); x[0] = xx[0];
	for (n=1, i=1; i < nn; i++) {
		if (!(xx[i] - xx[i-1] <= FLT_EPSILON)) {
			x[n] = xx[i]; f[n] = ff[i]; n++;
		}
	}
	spl->n = n; if (n <= 3) sf_error("insufficient samples");
	spl->x = x;
	spl->f = f;
	spl->linear = linear;

	if (linear) {
		A = sf_floatalloc(n-1);
		B = sf_floatalloc(n-1);

		for (i=0; i < n-1; i++) {
			A[i] = f[i];
			B[i] = (f[i+1] - f[i]) / (x[i+1] - x[i]);
		}
		spl->a = A;    spl->b = B;
		spl->c = NULL; spl->d = NULL;
	} else {
		h = sf_floatalloc(n-1);	b = sf_floatalloc(n-1);
		u = sf_floatalloc(n-2);	v = sf_floatalloc(n-2);
		A = sf_floatalloc(n-1);	B = sf_floatalloc(n-1);
		C = sf_floatalloc(n-1);	D = sf_floatalloc(n-1);

		for (i=0; i < n-1; i++) {
			h[i] = tmp = x[i+1] - x[i];
			b[i] = (f[i+1] - f[i]) / tmp;
		}
		for (u[0]=2.*(h[1]+h[0]) , v[0]=6.*(b[1]-b[0]), i=1; i < n-2; i++) {
			tmp = u[i-1];
			u[i] = 2. * (h[i+1] + h[i]) - h[i]*h[i] / tmp;
			v[i] = 6. * (b[i+1] - b[i]) - h[i]*v[i-1] / tmp;
		}
		mysgtsv(n-2,1,&h[1],u,&h[1],v,n-2);

		/* natural boudary */
		A[0] = f[0];
		B[0] = h[0]*v[0]/(-6.) + (f[1]-f[0])/h[0];
		C[0] = 0;
		D[0] = v[0]/(6.*h[0]);
		for (i=1; i < n-1; i++) {
			A[i] = f[i];
			B[i] = h[i]*v[i]/(-6.) + h[i]*v[i-1]/(-3.) + (f[i+1]-f[i])/h[i];
			C[i] = .5 * v[i-1];
			D[i] = (v[i] - v[i-1]) / (6.*h[i]);
		}

		spl->a = A; spl->b = B;
		spl->c = C;	spl->d = D;
		free(h); free(b); free(u); free(v);
	}

	return spl;
}

void spline_eval(/*@out@*/ float *g, /* [m] */
				 const float *y,     /* [m] sorted */
				 int m,
				 const Spl *spl)
/*< evaluate >*/
{
	int i,j,n = spl->n;
	float t;
	const float *x = spl->x,*f = spl->f,
		*a = spl->a,*b = spl->b,
		*c = spl->c,*d = spl->d;

	for (j=0, i=0; i < m; i++) {
		/* x[j] <= y[i] <= x[j+1] */
		if      (y[i] < x[0  ]) g[i] = f[0  ];
		else if (y[i] > x[n-1]) g[i] = f[n-1];
		else {
			while ((x[j]-y[i]) * (x[j+1]-y[i]) > 0) j++;
			t = y[i] - x[j];
			if (spl->linear) {
				g[i] = a[j] + b[j] * t;
			} else {
				g[i] = a[j] + b[j] * t + c[j] * t*t + d[j] * t*t*t;
			}
		}
	}
}

void spline_free(Spl *spl)
/*< clean >*/
{
	free(spl->a); free(spl->b);
	if (!spl->linear) {
		free(spl->c); free(spl->d);
	}
	free(spl);
}
