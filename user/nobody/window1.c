/* Splitting data into windows, 1-D */
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

#include "window1.h"

static float h,dw;

void window1_init (float h_in  /* window overlap */, 
		   float dw_in /* window size */)
/*< initialize >*/ 
{
    h = h_in;
    dw = dw_in;
}

int window1_apply (int iw           /* window number */, 
		   int w            /* window size */, 
		   const float* dat /* input data */, 
		   bool left        /* taper on the left */, 
		   bool right       /* taper on the right */, 
		   float *win       /* output window */)
/*< extract window >*/
{
    int i, i0;
    float gain;

    i0 = 0.5+iw*dw;
    for (i=0; i < w; i++) {
	if (left && i < h) {
	    gain = sinf(0.5*SF_PI*i/h);
	    win[i] = gain*gain*dat[i0+i];
	} else if (right && i >= w-h) {
	    gain =  sinf(0.5*SF_PI*(w-i-1)/h);
	    win[i] = gain*gain*dat[i0+i];
	} else {
	    win[i] = dat[i0+i];
	}
    }

    return i0;
}

/* 	$Id$	 */
