/* Chebyshev interpolation */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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
static float m, h, *c, *ci=NULL;

void chebyshev_init(int n1,                /* data size */
		    float xmin, float xmax /* axis range */)
/*< initialize >*/
{
    n = n1;

    c = sf_floatalloc(n);
    sf_cosft_init(n);

    m = (xmin+xmax)*0.5f;
    h = 2.0f/(xmax-xmin);
}

void cchebyshev_init(int n1,               /* data size */
		    float xmin, float xmax /* axis range */)
/*< initialize >*/
{
    n = n1;

    c = sf_floatalloc(n);
    ci = sf_floatalloc(n);
    sf_cosft_init(n);

    m = (xmin+xmax)*0.5f;
    h = 2.0f/(xmax-xmin);
}


void chebyshev_close(void)
/*< free allocated storage >*/
{
    free(c);
    sf_cosft_close();
}

void cchebyshev_close(void)
/*< free allocated storage >*/
{
    free(c);
    free(ci);
    sf_cosft_close();
}

void chebyshev_set(const float *d /* [n] data at Chebyshev points */)
/*< compute Chebyshev coefficients >*/
{
    int i;

    for (i=0; i < n; i++) {
	c[i] = d[i];
    }
    sf_cosft_inv(c,0,1);
}

void cchebyshev_set(const sf_complex *d /* [n] data at Chebyshev points */)
/*< compute Chebyshev coefficients >*/
{
    int i;

    for (i=0; i < n; i++) {
	c[i] = crealf(d[i]);
	ci[i] = cimagf(d[i]);
    }
    sf_cosft_inv(c,0,1);
    sf_cosft_inv(ci,0,1);
}

void chebyshev_poly(float *c2)
/*< return chebyshev coefficients >*/
{
    int i;

    for (i=0; i < n; i++) {
	c2[i] = c[i];
    }
}

void cchebyshev_poly(sf_complex *c2)
/*< return chebyshev coefficients >*/
{
    int i;

    for (i=0; i < n; i++) {
	c2[i] = sf_cmplx(c[i],ci[i]);
    }
}

float chebyshev(float x)
/*< interpolate >*/
{
    int i;
    float c0, c1, c2;

    x = (x-m)*h;
    c1 = 0.0f;
    c0 = c[n-1];
    for (i=n-2; i > 0; i--) {
	c2 = c1;
	c1 = c0;
	c0 = 2*(c[i] + x*c0) - c2;
    }
    c0 = c[0] + x*c0 - c1;
    return c0;
}

sf_complex cchebyshev(float x)
/*< interpolate >*/
{
    int i;
    float cr, cj, c1, c2;

    x = (x-m)*h;
    
    c1 = 0.0f;
    cr = c[n-1];
    for (i=n-2; i > 0; i--) {
	c2 = c1;
	c1 = cr;
	cr = 2*(c[i] + x*cr) - c2;
    }
    cr = c[0] + x*cr - c1;

    c1 = 0.0f;
    cj = ci[n-1];
    for (i=n-2; i > 0; i--) {
	c2 = c1;
	c1 = cj;
	cj = 2*(ci[i] + x*cj) - c2;
    }
    cj = ci[0] + x*cj - c1;
    
    return sf_cmplx(cr,cj);
}

