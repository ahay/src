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

#include <rsf.h>

#include "spline.h"
#ifndef _spline_h

typedef struct{
	int n;
	float *x,*f;       /*[n]   data f(x)    */
	float *a,*b,*c,*d; /*[n-1] coefficients */
} Spl;
/*^*/

#endif

void trid(/*@out@*/ float *x, /* [n]   */
		  const float *d,     /* [n]   */
		  const float *u,     /* [n-1] */
		  const float *l,     /* [n-1] */
		  const float *y,     /* [n]   */
		  int n)
/* tridiagonal system solver */
{
	int i;
	float *dd,*yy,tmp;

	dd = sf_floatalloc(n);
	yy = sf_floatalloc(n);

	/* naive elimination */
	for (dd[0]=d[0], yy[0]=y[0], i=1; i < n; i++) {
		tmp = l[i-1] / dd[i-1];
		dd[i] = d[i] - tmp * u[i-1];
		yy[i] = y[i] - tmp *yy[i-1];
	}
	/* back substitution */
	for (x[n-1]=yy[n-1]/dd[n-1], i=n-2; i >= 0; i--) {
		x[i] =(yy[i] - u[i]*x[i+1]) / dd[i];
	}
}

Spl *spline_init(const float *f, /* [n] */
				 const float *x, /* [n] */
				 int n)
/*< initialize >*/
{
	int i;
	float tmp,*h,*b,*u,*v,*z,*A,*B,*C,*D;
	Spl *spl;

	h = sf_floatalloc(n - 1);	b = sf_floatalloc(n - 1);
	u = sf_floatalloc(n - 2);	v = sf_floatalloc(n - 2);
	z = sf_floatalloc(n);
	A = sf_floatalloc(n - 1);	B = sf_floatalloc(n - 1);
	C = sf_floatalloc(n - 1);	D = sf_floatalloc(n - 1);

	for (i=0; i < n-1; i++) {
		h[i] = tmp = x[i+1] - x[i];
		b[i] = (f[i+1] - f[i]) / tmp;
	}
	for (u[0]=2.*(h[1]+h[0]) , v[0]=6.*(b[1]-b[0]), i=1; i < n-2; i++) {
		u[i] = 2. * (h[i+1] + h[i]) - h[i]*h[i]/u[i-1];
		v[i] = 6. * (b[i+1] - b[i]) - h[i]*v[i-1]/u[i-1];
	}
	trid(&z[1],u,&h[1],&h[1],v,n-2);
	/* natural boudary */
	z[0] = 0.; z[n-1] = 0;
#ifdef DEBUG
	for (i=0; i < n; i++) sf_warning("z[%d] = %g",i,z[i]);
#endif

	for (i=0; i < n-1; i++) {
		A[i] = f[i];
		B[i] = h[i]*z[i+1]/(-6.) + h[i]*z[i]/(-3.) + (f[i+1]-f[i])/h[i];
		C[i] = .5 * z[i];
		D[i] = (z[i+1] - z[i]) / (6.*h[i]);
#ifdef DEBUG
		sf_warning("A[%d]=%g,B[%d]=%g,C[%d]=%g,D[%d]=%g",i,A[i],i,B[i],i,C[i],i,D[i]);
#endif
	}

	spl = sf_alloc(1,sizeof(*spl));
	spl->n = n;
	spl->x = sf_floatalloc(n); spl->f = sf_floatalloc(n);
	for (i=0; i<n; i++) {spl->x[i] = x[i]; spl->f[i] = f[i];}
	spl->a = A; spl->b = B; spl->c = C; spl->d = D;
	return spl;
}

void spline_eval(/*@out@*/ float *g, /* [m] g(y) */
				 const float *y,     /* [m] y, y[i] <= y[i+1] */
				 int m,
				 const Spl *spl)
/*< evaluate >*/
{
	int i,j,n;
	float t,tt,ttt,*x,*f,*a,*b,*c,*d;
	
	x = spl->x; f = spl->f; n = spl->n;
	a = spl->a; b = spl->b; c = spl->c; d = spl->d;

	/* x[j] <= y[i] <= x[j+1] */
	for (j=0, i=0; i < m; i++) {
		if      (y[i] < x[0  ]) g[i] = f[0  ];
		else if (y[i] > x[n-1]) g[i] = f[n-1];
		else {
			while ((x[j]-y[i]) * (x[j+1]-y[i]) > 0) j++;
			t = y[i] - x[j]; tt = t*t; ttt = tt*t;
			g[i] = a[j] + b[j]*t + c[j]*tt + d[j]*ttt;
		}
	}
}

void spline_free(Spl *spl)
/*< clean >*/
{
	free(spl->f); free(spl->x);
	free(spl->a); free(spl->b); free(spl->c); free(spl->d);
	free(spl);
}
