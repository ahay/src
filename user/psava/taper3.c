/* Tapering (2d and 3d) */
/*
  Copyright (C) 2007 Colorado School of Mines
   
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

#include "taper3.h"

#include "weutil.h"
/*^*/

/*------------------------------------------------------------*/
tap3d taper_init(int n1_,
		 int n2_,
		 int n3_,
		 int   nt1_, 
		 int   nt2_, 
		 int   nt3_ /* taper lengths */,
		 bool  b1_,
		 bool  b2_,
		 bool  b3_)
/*< 3-D initialize >*/
{
    int it;
    float gain;

    tap3d tap;
    tap = (tap3d) sf_alloc(1,sizeof(*tap));

    tap->n1 =n1_;
    tap->n2 =n2_;
    tap->n3 =n3_;

    tap->nt1=nt1_;
    tap->nt2=nt2_;
    tap->nt3=nt3_;

    tap->b1 =b1_;
    tap->b2 =b2_;
    tap->b3 =b3_;

    if (tap->nt1 > 0) {
	tap->tap1 = sf_floatalloc(tap->nt1);
	for (it=0; it < tap->nt1; it++) {
	    gain = sinf(0.5*SF_PI*it/tap->nt1);
	    tap->tap1[it]=(1+gain)/2.;
	}
    }
    if (tap->nt2 > 0) {
	tap->tap2 = sf_floatalloc(tap->nt2);
	for (it=0; it < tap->nt2; it++) {
	    gain = sinf(0.5*SF_PI*it/tap->nt2);
	    tap->tap2[it]=(1+gain)/2.;
	}
    }
    if (tap->nt3 > 0) {
	tap->tap3 = sf_floatalloc(tap->nt3);
	for (it=0; it < tap->nt3; it++) {
	    gain = sinf(0.5*SF_PI*it/tap->nt3);
	    tap->tap3[it]=(1+gain)/2.;
	}
    }

    return tap;
}

/*------------------------------------------------------------*/
void taper_report(tap3d tap)
/*< report taper parameters >*/
{

    sf_warning("TAPER REPORT");
    sf_warning(" n1=%d  n2=%d  n3=%d ",tap->n1, tap->n2, tap->n3);
    sf_warning("nt1=%d nt2=%d nt3=%d",tap->nt1,tap->nt2,tap->nt3);
    sf_warning(" b1=%d  b2=%d  b3=%b ",tap->b1, tap->b2, tap->b3);

}

/*------------------------------------------------------------*/
void taper2d_close(tap3d tap)
/*< 2-D free allocated storage >*/
{
    if (tap->nt1 > 0) free(tap->tap1);
    if (tap->nt2 > 0) free(tap->tap2);
}

/*------------------------------------------------------------*/
void taper3d_close(tap3d tap)
/*< 3-D free allocated storage >*/
{
    if (tap->nt1 > 0) free(tap->tap1);
    if (tap->nt2 > 0) free(tap->tap2);
    if (tap->nt3 > 0) free(tap->tap3);
}

/*------------------------------------------------------------*/
void taper2d(sf_complex **tt  /* [n2][n1] tapered array (in and out) */,
	     tap3d tap)
/*< 2-D taper >*/
{
    int it,i2,i1;
    float gain;

    for (it=0; it < tap->nt2; it++) {
	gain = tap->tap2[it];
	for (i1=0; i1 < tap->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    if (tap->b2) tt[        it  ][i1] *= gain;
	    ;            tt[tap->n2-it-1][i1] *= gain;
#else
	    if (tap->b2) tt[        it  ][i1] = sf_crmul(tt[        it  ][i1],gain);
	    ;            tt[tap->n2-it-1][i1] = sf_crmul(tt[tap->n2-it-1][i1],gain);
#endif
	}
    }

    for (it=0; it < tap->nt1; it++) {
	gain = tap->tap1[it];
	for (i2=0; i2 < tap->n2; i2++) {
#ifdef SF_HAS_COMPLEX_H
	    if (tap->b1) tt[i2][        it  ] *= gain;
	    ;            tt[i2][tap->n1-it-1] *= gain;
#else
	    if (tap->b1) tt[i2][        it  ] = sf_crmul(tt[i2][        it  ],gain);
	    ;            tt[i2][tap->n1-it-1] = sf_crmul(tt[i2][tap->n1-it-1],gain);
#endif
	}
    }
}

/*------------------------------------------------------------*/
void taper3d(sf_complex*** tt /* [n3][n2][n1] tapered array (inout) */,
	     tap3d tap)
/*< 3-D taper >*/
{
    int it,i1,i2,i3;
    float gain;

    for (it=0; it < tap->nt3; it++) {
	gain = tap->tap3[it];
	for     (i2=0; i2 < tap->n2; i2++) {
	    for (i1=0; i1 < tap->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		if(tap->b3) tt[        it  ][i2][i1] *= gain;
		;           tt[tap->n3-it-1][i2][i1] *= gain;
#else
		if(tap->b3)  tt[        it  ][i2][i1] = 
		    sf_crmul(tt[        it  ][i2][i1],gain);
		;            tt[tap->n3-it-1][i2][i1] = 
		    sf_crmul(tt[tap->n3-it-1][i2][i1],gain);
#endif
	    }
	}
    }

    for (it=0; it < tap->nt2; it++) {
	gain = tap->tap2[it];
	for     (i3=0; i3 < tap->n3; i3++) {
	    for (i1=0; i1 < tap->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		if(tap->b2) tt[i3][        it  ][i1] *= gain;
		;           tt[i3][tap->n2-it-1][i1] *= gain;
#else
		if(tap->b2) tt[i3][        it  ][i1] = 
		   sf_crmul(tt[i3][        it  ][i1],gain);
		;           tt[i3][tap->n2-it-1][i1] = 
		   sf_crmul(tt[i3][tap->n2-it-1][i1],gain);
#endif
	    }
	}
    }

    for (it=0; it < tap->nt1; it++) {
	gain = tap->tap1[it];
	for     (i3=0; i3 < tap->n3; i3++) {
	    for (i2=0; i2 < tap->n2; i2++) {
#ifdef SF_HAS_COMPLEX_H
		if(tap->b1) tt[i3][i2][        it  ] *= gain;
		;           tt[i3][i2][tap->n1-it-1] *= gain;
#else
		if(tap->b1) tt[i3][i2][        it  ] = 
		   sf_crmul(tt[i3][i2][        it  ],gain);
		;           tt[i3][i2][tap->n1-it-1] = 
		   sf_crmul(tt[i3][i2][tap->n1-it-1],gain);
#endif
	    }
	}
    }
}
