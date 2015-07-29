/* Implicit step in Chebyshev-tau method. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

static int n;
static float *e, *f, *u, *a, *c, *g;

void tricheb_init (int nd /* data size */)
/*< initialize >*/
{
    n = nd; 
    e = sf_floatalloc(n);
    f = sf_floatalloc(n);
    u = sf_floatalloc(n);
    a = sf_floatalloc(n);
    c = sf_floatalloc(n);
    g = sf_floatalloc(n);
    a[n-1] = 0.; 
    g[n-1] = 1.;
}

void tricheb_close(void)
/*< free allocated storage >*/
{
    free(e);
    free(f);
    free(u);
    free(a);
    free(c);
    free(g);
}

void tricheb_step (float s, float z)
/*< one step >*/
{
    int i;
    float t;

    c[n-1] = 0.5*z/n;
    e[n-1] = -c[n-1];
    t = z; 
    u[n-1] = z;
    for (i=n-2; i >= 0; i--) {
	a[i] = -0.5*z/(i+1);
	c[i] = +0.5*z/(i+1);
	if (s < 0.) t = -t;
	g[i]   = 1.+a[i]*e[i+1];   
	u[i] = t +u[i+1]*e[i+1];
	e[i] = -c[i]/g[i];
    }
    c[0] = z;
}

void tricheb_apply (float *d, const float *q)
/*< forward >*/
{
    int i;

    for (i=0; i < n-1; i++) {
	d[i] = q[i+1] - a[i]*q[i+2] - c[i]*q[i];
    }
    d[n-1] = q[n] - c[n-1]*q[n-1];
}

void tricheb_solve (const float *d, float *q)
/*< inverse >*/
{
    int i;
    float dp;

    f[n-1] = d[n-1];
    for (i=n-2; i >=0; i--) {
	f[i] = (d[i]-a[i]*f[i+1])/g[i];
    }

    dp = 0.;
    for (i=0; i < n; i++) {
	dp += u[i]*f[i];
    }
    q[0] = f[0]+dp/(g[0]-u[0]);
    for (i=1; i < n; i++) {
	q[i] = e[i] * q[i-1] + f[i];
    }
}

