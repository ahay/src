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
		    const float *d,        /* [n] data at Chebyshev points */
		    float xmin, float xmax /* axis range */)
/*< initialize >*/
{
    int i;

    n = n1;

    c = sf_floatalloc(n);
    for (i=0; i < n; i++) {
	c[i] = d[i];
    }

    sf_cosft_init(n);
    sf_cosft_inv(c,0,1);
    sf_cosft_close();

    for (i=1; i < n-1; i++) {
	c[i] *= 2;
    }

    m = (xmin+xmax)*0.5f;
    h = 2.0f/(xmax-xmin);
}

float chebyshev(float x)
/*< interpolate >*/
{
    int i;

    float w0, w1, w2;

    x = (x-m)*h;

    w1=w2=0.0f;
    
    for (i=n-1; i >=0; i--) {
	w0 = c[i] + 2*x*w1 - w2;
	w2 = w1;
	w1 = w0;
    }
    return (w1-x*w2);
}



