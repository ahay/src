/* Hilbert transform FIR filter. */
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

#include "hilbert.h"

#include "alloc.h"

static float c, c2, *h;
static int n, nt;

void sf_hilbert_init(int nt1  /* trace length */, 
		     int n1   /* transform length */, 
		     float c1 /* filter parameter */)
/*< initialize >*/
{
    n = n1;
    nt = nt1;
    c = 1./(2*sqrtf(c1));
    c2 = c*c;
    h = sf_floatalloc(nt);
}

void sf_hilbert_close(void)
/*< free allocated storage >*/
{
    free(h);
}

void sf_hilbert (const float* trace, float* trace2)
/*< transform >*/
{
    int i, it;
    
    for (it=0; it < nt; it++) {
	h[it] = trace[it];
    }

    for (i=n; i >= 1; i--) {
	trace2[0] = h[0] + (h[2]-2*h[1]+h[0])*c2;
	trace2[1] = h[1] + (h[3]-3*h[1]+2*h[0])*c2;
	for (it=2; it < nt-2; it++) {
	    trace2[it] = h[it]+(h[it+2]-2.*h[it]+h[it-2])*c2;
	}
	trace2[nt-2] = h[nt-2] + (h[nt-4]-3*h[nt-2]+2*h[nt-1])*c2;
	trace2[nt-1] = h[nt-1] + (h[nt-3]-2*h[nt-2]+h[nt-1])*c2;

	for (it=0; it < nt; it++) {
	    h[it] = trace[it] + trace2[it]*(2*i-1)/(2*i);
	}
    }

    trace2[0] = 2.*(h[0]-h[1])*c;
    for (it=1; it < nt-1; it++) {
	trace2[it] = (h[it-1]-h[it+1])*c;
    }
    trace2[nt-1] = 2.*(h[nt-2]-h[nt-1])*c;
}

void sf_hilbert4 (const float* trace, float* trace2)
/*< transform - kind 4 filter >*/
{
    int i, it;
    
    for (it=0; it < nt; it++) {
	h[it] = trace[it];
    }

    for (i=n; i >= 1; i--) {
	for (it=1; it < nt-1; it++) {
	    trace2[it] = h[it]+(h[it+1]-2.*h[it]+h[it-1])*c2;
	}
	trace2[0] = trace2[1];
	trace2[nt-1] = trace2[nt-2];

	for (it=0; it < nt; it++) {
	    h[it] = trace[it] + trace2[it]*(2*i-1)/(2*i);
	}
    }

    for (it=0; it < nt-1; it++) {
	trace2[it] = (h[it]-h[it+1])*c;
    }
    trace2[nt-1] = trace2[nt-2];
}

