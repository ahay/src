/* Extract window from data, 3-D */
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "window3.h"

static int *w, *nw, *n;
static float *h;

void window3_init (int* w_in   /* window size [3] */, 
		   int* nw_in  /* number of windows [3] */, 
		   int* n_in   /* data size [3] */, 
		   float* h_in /* overlap [3] */)
/*< initialize >*/
{
    w = w_in; 
    nw = nw_in; 
    n = n_in; 
    h = h_in;
}

void window3_apply (const int* i           /* window number [3] */, 
		    float*** dat           /* input data */, 
		    bool near, bool far, 
		    bool left, bool right, 
		    bool top, bool bottom  /* tapering */,
		    int* i0                /* window start [3] */, 
		    float*** win           /* output window */)
/*< extract window >*/
{
    int i1, i2, i3, j;
    float gain1, gain2, gain3;

    for (j=0; j < 3; j++) {
	i0[j] = 0.5+(n[j]-w[j])*i[j]/(nw[j]-1.);
    }

    for (i3=0; i3 < w[2]; i3++) {
	if (near && i3+1 <= h[2]) {
	    gain3 = sinf(0.5*SF_PI*i3/h[2]);
	} else if (far && i3 >= w[2]-h[2]) {
	    gain3 = sinf(0.5*SF_PI*(w[2]-i3-1)/h[2]);
	} else {
	    gain3 = 1.;
	}
	for (i2=0; i2 < w[1]; i2++) {
	    if (left && i2+1 <= h[1]) {
		gain2 = sinf(0.5*SF_PI*i2/h[1]);
	    } else if (right && i2 >= w[1]-h[1]) {
		gain2 = sinf(0.5*SF_PI*(w[1]-i2-1)/h[1]);
	    } else {
		gain2 = 1.;
	    }
	    
	    gain2 *= gain3;

	    for (i1=0; i1 < w[0]; i1++) {
		if (top && i1+1 <= h[0]) {
		    gain1 = sinf(0.5*SF_PI*i1/h[0]);
		} else if (bottom && i1 >= w[0]-h[0]) {
		    gain1 = sinf(0.5*SF_PI*(w[0]-i1-1)/h[0]);
		} else {
		    gain1 = 1.;
		}
		
		gain1 *= gain2;
		gain1 *= gain1;
          
		win[i3][i2][i1] = gain1*dat[i0[2]+i3][i0[1]+i2][i0[0]+i1];
	    }
	}
    }
}

/* 	$Id$	 */
