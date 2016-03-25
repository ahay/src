/* First derivative FIR filter. */
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

/* From "Closed-form design of maximally flat FIR Hilbert transformers,
 * differentiators, and fractional delayers by power series expansion" by
 * Soo-Chang Pei and Peng-Hua Wang, IEEE Trans. on Circuits and Systems - Part
 * I: Fundamental theory and applications, v. 48, No. 4, 2001, 389-398. */

#include <math.h>

#include "deriv.h"
#include "alloc.h"

static float c, c2, *h;
static int n, nt;

void sf_deriv_init(int nt1  /* transform length */, 
		   int n1   /* trace length */, 
		   float c1 /* filter parameter */)
/*< initialize >*/
{
    n = n1;
    nt = nt1;
    c = 1./(2*sqrtf(c1));
    c2 = c*c; c2=c2*1.0;
    h = sf_floatalloc(nt);
}

void sf_deriv_close(void)
/*< free allocated storage >*/
{
    free(h);
}

void sf_deriv (const float* trace, float* trace2)
/*< derivative operator >*/
{
    int i, it;
    
    for (it=0; it < nt; it++) {
	h[it] = trace[it];
    }

    for (i=n; i >= 1; i--) {
	for (it=1; it < nt-1; it++) {
	    trace2[it] = h[it]-0.5*(h[it+1]+h[it-1]);
	}
	trace2[0] = trace2[1];
	trace2[nt-1] = trace2[nt-2];

	for (it=0; it < nt; it++) {
	    h[it] = trace[it] + trace2[it]*i/(2*i+1);
	}
    }

    trace2[0] = h[1]-h[0];
    for (it=1; it < nt-1; it++) {
	trace2[it] = 0.5*(h[it+1]-h[it-1]);
    }
    trace2[nt-1] = h[nt-1]-h[nt-2];
}

