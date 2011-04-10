/* Extract window from data, 2-D */
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

#include "window2.h"

static int *w, *nw, *n;
static float *h;

void window2_init (int* w_in   /* window size [2] */, 
		   int* nw_in  /* number of windows [2] */, 
		   int* n_in   /* data size [2] */, 
		   float* h_in /* window overlap [2] */)
/*< initialize >*/
{
    w = w_in; 
    nw = nw_in; 
    n = n_in; 
    h = h_in;
}

void window2_apply (const int* i /* window number [2] */, 
		    float** dat  /* input data */, 
		    bool left    /* left taper */, 
		    bool right   /* right taper */, 
		    bool top     /* top taper */, 
		    bool bottom  /* bottom taper */,
		    int* i0      /* window start [2] */, 
		    float** win  /* output window */)
/*< extract window >*/
{
    int i1, i2, j;
    float gain1, gain2;

    for (j=0; j < 2; j++) {
	i0[j] = 0.5+(n[j]-w[j])*i[j]/(nw[j]-1.);
    }

    for (i2=0; i2 < w[1]; i2++) {
	if (left && i2+1 <= h[1]) {
	    gain2 = sinf(0.5*SF_PI*i2/h[1]);
	} else if (right && i2 >= w[1]-h[1]) {
	    gain2 = sinf(0.5*SF_PI*(w[1]-i2-1)/h[1]);
	} else {
	    gain2 = 1.;
	}
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
          
	    win[i2][i1] = gain1*dat[i0[1]+i2][i0[0]+i1];
	}
    }
}

/* 	$Id$	 */
