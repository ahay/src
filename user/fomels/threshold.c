/* sharpening */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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
static float *ww, s, dmin, dmax;

void threshold_init(int n1, float scale) 
/*< initialize >*/
{
    n = n1;
    s = scale;
    
    ww = sf_floatalloc(n);
}

void threshold_set(const float *data)
/*< find threshold >*/
{
    int i;
    float q1, q3;

    for (i=0; i < n; i++) {
	ww[i] = data[i];
    }

    q1 = sf_quantile(ceilf(0.25*n),n,ww);
    q3 = sf_quantile(floorf(0.75*n),n,ww);

    dmin = q1 - s*(q3-q1);
    dmax = q3 + s*(q3-q1);
}

void threshold_close(void)
/*< free allocated storage >*/
{
    free(ww);
}

void threshold(float *data)
/*< apply thresholding >*/
{
    int i;
    float d;

    for (i=0; i < n; i++) {
	d = data[i];
	if (d < dmin) {
	    data[i] -= dmin;
	} else if (d > dmax) {
	    data[i] -= dmax;
	} else {
	    data[i] = 0.0f;
	}
    }
}

