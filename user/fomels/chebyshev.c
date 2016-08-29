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
static float m, h, *c;

void chebyshev_init(int n1,                 /* data size */
		    float xmin, float xmax /* axis range */)
/*< initialize >*/
{
    n = n1;

    c = sf_floatalloc(n);

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

void chebyshev_set(const float *d /* [n] data at Chebyshev points */)
/*< compute Chebyshev coefficients >*/
{
    int i;

    for (i=0; i < n; i++) {
	c[i] = d[i];
    }
    sf_cosft_inv(c,0,1);
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




