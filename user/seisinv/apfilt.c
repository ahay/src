/* All-pass plane-wave destruction filter coefficients */
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

#include <rsf.h>

#include "apfilt.h"

static int n;
static float *b;

void apfilt_init(int nw /* filter order */)
/*< initialize >*/
{
    int j, k;
    double bk;

    n = nw*2;
    b = sf_floatalloc(n+1);

    for (k=0; k <= n; k++) {
	bk = 1.0;
	for (j=0; j < n; j++) {
	    if (j < n-k) {
		bk *= (k+j+1.0)/(2*(2*j+1)*(j+1));
	    } else {
		bk *= 1.0/(2*(2*j+1));
	    }
	}
	b[k] = bk;
    }
}

void apfilt_close(void)
/*< free allocated storage >*/
{
    free(b);
}

void passfilter (float p  /* slope */, 
		 float* a /* output filter [n+1] */)
/*< find filter coefficients >*/
{
    int j, k;
    float ak;
    
    for (k=0; k <= n; k++) {
	ak = b[k];
	for (j=0; j < n; j++) {
	    if (j < n-k) {
		ak *= (n-j-p);
	    } else {
		ak *= (p+j+1);
	    }
	}
	a[k] = ak;
    }
}

void aderfilter (float p  /* slope */, 
		 float* a /* output filter [n+1] */)
/*< find coefficients for filter derivative >*/
{

    int i, j, k;
    float ak, ai;
    
    for (k=0; k <= n; k++) {
	ak = 0.;
	for (i=0; i < n; i++) {
	    ai = -1.0;
	    for (j=0; j < n; j++) {
		if (j != i) {			
		    if (j < n-k) {
			ai *= (n-j-p);
		    } else {
			ai *= (p+j+1);
		    }
		} else if (j < n-k) {
		    ai *= (-1);
		}
	    }
	    ak += ai;
	}
	a[k] = ak*b[k];
    }
}

/* 	$Id: apfilt.c 4781 2009-09-25 04:01:38Z sfomel $	 */
