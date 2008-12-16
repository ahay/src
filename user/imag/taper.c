/* Tapering */
/*
  Copyright (C) 2006 Colorado School of Mines
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
static int  n1,  n2,  n3;
static bool b1,  b2,  b3;
static float *tap1=NULL, *tap2=NULL, *tap3=NULL;

void taper2_init(int n2_, int n1_ /* cube dimensions */,
		 int m2_, int m1_ /* taper lengths */,
		 bool b2_, bool b1_)
/*< 2-D initialize >*/
{

    int it;
    float gain;

    nt1=m1_;
    nt2=m2_;
    
    n1 =n1_;
    n2 =n2_;

    b1 =b1_;
    b2 =b2_;

    if (nt1 > 0) {
	tap1 = sf_floatalloc(nt1);
	for (it=0; it < nt1; it++) {
	    gain = sinf(0.5*SF_PI*it/nt1);
	    tap1[it]=(1+gain)/2.;
	}
    }
    if (nt2 > 0) {
	tap2 = sf_floatalloc(nt2);
	for (it=0; it < nt2; it++) {
	    gain = sinf(0.5*SF_PI*it/nt2);
	    tap2[it]=(1+gain)/2.;
	}
    }
}

void taper2_close(void)
/*< 2-D free allocated storage >*/
{
    free(tap1);
    free(tap2);
}

void taper3_init(
    int n3_, 
    int n2_, 
    int n1_ /* cube dimensions */,
    int m3_, 
    int m2_, 
    int m1_ /* taper lengths */,
    bool b3_,
    bool b2_,
    bool b1_)
/*< 3-D initialize >*/
{
    int it;
    float gain;

    nt1=m1_;
    nt2=m2_;
    nt3=m3_;

    n1 =n1_;
    n2 =n2_;
    n3 =n3_;

    b1 =b1_;
    b2 =b2_;
    b3 =b3_;

    if (nt1 > 0) {
	tap1 = sf_floatalloc(nt1);
	for (it=0; it < nt1; it++) {
	    gain = sinf(0.5*SF_PI*it/nt1);
	    tap1[it]=(1+gain)/2.;
	}
    }
    if (nt2 > 0) {
	tap2 = sf_floatalloc(nt2);
	for (it=0; it < nt2; it++) {
	    gain = sinf(0.5*SF_PI*it/nt2);
	    tap2[it]=(1+gain)/2.;
	}
    }
    if (nt3 > 0) {
	tap3 = sf_floatalloc(nt3);
	for (it=0; it < nt3; it++) {
	    gain = sinf(0.5*SF_PI*it/nt3);
	    tap3[it]=(1+gain)/2.;
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

void taper2(sf_complex** tt  /* [n2][n1] tapered array (in and out) */)
/*< 2-D taper >*/
{
    int it,i2,i1;
    float gain;

    for (it=0; it < nt2; it++) {
	gain = tap2[it];
	for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    if (b2) tt[   it  ][i1] *= gain;
	    ;       tt[n2-it-1][i1] *= gain;
#else
	    if (b2) tt[   it  ][i1] = sf_crmul(tt[   it  ][i1],gain);
	    ;       tt[n2-it-1][i1] = sf_crmul(tt[n2-it-1][i1],gain);
#endif
	}
    }

    for (it=0; it < nt1; it++) {
	gain = tap1[it];
	for (i2=0; i2 < n2; i2++) {
#ifdef SF_HAS_COMPLEX_H
	    if (b1) tt[i2][   it  ] *= gain;
	    ;       tt[i2][n1-it-1] *= gain;
#else
	    if (b1) tt[i2][   it  ] = sf_crmul(tt[i2][   it  ],gain);
	    ;       tt[i2][n1-it-1] = sf_crmul(tt[i2][n1-it-1],gain);
#endif
	}
    }
}

void taper3(sf_complex*** tt /* [n3][n2][n1] tapered array (inout) */)
/*< 3-D taper >*/
{
    int it,i1,i2,i3;
    float gain;

    for (it=0; it < nt3; it++) {
	gain = tap3[it];
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		if(b3) tt[   it  ][i2][i1] *= gain;
		;      tt[n3-it-1][i2][i1] *= gain;
#else
		if(b3) tt[   it  ][i2][i1] = 
			   sf_crmul(tt[   it  ][i2][i1],gain);
		;      tt[n3-it-1][i2][i1] = 
			   sf_crmul(tt[n3-it-1][i2][i1],gain);
#endif
	    }
	}
    }

    for (it=0; it < nt2; it++) {
	gain = tap2[it];
	for (i3=0; i3 < n3; i3++) {
	    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		if(b2) tt[i3][   it  ][i1] *= gain;
		;      tt[i3][n2-it-1][i1] *= gain;
#else
		if(b2) tt[i3][   it  ][i1] = 
			   sf_crmul(tt[i3][   it  ][i1],gain);
		;      tt[i3][n2-it-1][i1] = 
			   sf_crmul(tt[i3][n2-it-1][i1],gain);
#endif
	    }
	}
    }

    for (it=0; it < nt1; it++) {
	gain = tap1[it];
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
#ifdef SF_HAS_COMPLEX_H
		if(b1) tt[i3][i2][   it  ] *= gain;
		;      tt[i3][i2][n1-it-1] *= gain;
#else
		if(b1) tt[i3][i2][   it  ] = 
			   sf_crmul(tt[i3][i2][   it  ],gain);
		;      tt[i3][i2][n1-it-1] = 
			   sf_crmul(tt[i3][i2][n1-it-1],gain);
#endif
	    }
	}
    }
}
