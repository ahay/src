/* Tapering */
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

#include "taper.h"

static int nt1, nt2, nt3;
static float *tap1=NULL, *tap2=NULL, *tap3=NULL;

void taper2_init(int n2, int n1 /* taper lengths */)
/*< 2-D initialize >*/
{
    int it;
    float gain;

    nt1=n1;
    nt2=n2;
    if (nt1 > 0) {
	tap1 = sf_floatalloc(nt1);
	for (it=0; it < nt1; it++) {
	    gain = sinf(0.5*SF_PI*it/nt1);
	    tap1[it]=gain*gain;
	}
    }
    if (nt2 > 0) {
	tap2 = sf_floatalloc(nt2);
	for (it=0; it < nt2; it++) {
	    gain = sinf(0.5*SF_PI*it/nt2);
	    tap2[it]=gain*gain;
	}
    }
}

void taper2_close(void)
/*< 2-D free allocated storage >*/
{
    free(tap1);
    free(tap2);
}

void taper3_init(int n3, int n2, int n1 /* taper lengths */)
/*< 3-D initialize >*/
{
    int it;
    float gain;

    nt1=n1;
    nt2=n2;
    nt3=n3;
    if (nt1 > 0) {
	tap1 = sf_floatalloc(nt1);
	for (it=0; it < nt1; it++) {
	    gain = sinf(0.5*SF_PI*it/nt1);
	    tap1[it]=gain*gain;
	}
    }
    if (nt2 > 0) {
	tap2 = sf_floatalloc(nt2);
	for (it=0; it < nt2; it++) {
	    gain = sinf(0.5*SF_PI*it/nt2);
	    tap2[it]=gain*gain;
	}
    }
    if (nt3 > 0) {
	tap3 = sf_floatalloc(nt3);
	for (it=0; it < nt3; it++) {
	    gain = sinf(0.5*SF_PI*it/nt3);
	    tap2[it]=gain*gain;
	}
    }
}

void taper3_close(void)
/*< 2-D free allocated storage >*/
{
    if (nt1 > 0) free(tap1);
    if (nt2 > 0) free(tap2);
    if (nt3 > 0) free(tap3);
}

void taper2(bool beg2, bool beg1  /* taper in the beginning  */, 
	    int n2, int n1        /* dimensions */, 
	    float complex** tt    /* [n2][n1] tapered array (in and out) */)
/*< 2-D taper >*/
{
    int it,i2,i1;
    float gain;

    for (it=0; it < nt2; it++) {
	gain = tap2[it];
	for (i1=0; i1 < n1; i1++) {
	    if (beg2) tt[it][i1] *= gain;
	    tt[n2-it-1][i1] *= gain;
	}
    }

    for (it=0; it < nt1; it++) {
	gain = tap1[it];
	for (i2=0; i2 < n2; i2++) {
	    if (beg1) tt[i2][it] *= gain;
	    tt[i2][n1-it-1] *= gain;
	}
    }
}

void taper3(bool b3, bool b2, bool b1 /* taper in the beginning  */, 
	    int n3, int n2, int n1    /* dimensions */, 
	    float complex*** tt       /* [n3][n2][n1] tapered array (inout) */)
/*< 3-D taper >*/
{
    int it,i1,i2,i3;
    float gain;

    for (it=0; it < nt3; it++) {
	gain = tap3[it];
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		if(b3) tt[   it  ][i2][i1] *= gain;
		;      tt[n3-it-1][i2][i1] *= gain;
	    }
	}
    }

    for (it=0; it < nt2; it++) {
	gain = tap2[it];
	for (i3=0; i3 < n3; i3++) {
	    for (i1=0; i1 < n1; i1++) {
		if(b2) tt[i3][   it  ][i1] *= gain;
		;      tt[i3][n2-it-1][i1] *= gain;
	    }
	}
    }

    for (it=0; it < nt1; it++) {
	gain = tap1[it];
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		if(b1) tt[i3][i2][   it  ] *= gain;
		;      tt[i3][i2][n1-it-1] *= gain;
	    }
	}
    }
}
