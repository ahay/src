/* Tapering (2d and 3d) */

/*
  Copyright (C) 2010 Colorado School of Mines
   
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

#include "weiutl.h"
/*^*/

#include "weitap.h"

/*------------------------------------------------------------*/
weitap3d weitap_init(weicub3d cub)
/*< taper initialize >*/
{
    int it;
    float gain;

    weitap3d tap;
    tap = (weitap3d) sf_alloc(1,sizeof(*tap));

    tap->n1 =sf_n(cub->amx);
    tap->n2 =sf_n(cub->amy);

    tap->nt1=SF_MIN(cub->tmx,sf_n(cub->amx)-1);
    tap->nt2=SF_MIN(cub->tmy,sf_n(cub->amy)-1);

    tap->b1 =true;
    tap->b2 =true;

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

    return tap;
}

/*------------------------------------------------------------*/
void weitap_close(weitap3d tap)
/*< 2-D free allocated storage >*/
{
    if (tap->nt1 > 0) free(tap->tap1);
    if (tap->nt2 > 0) free(tap->tap2);
}

/*------------------------------------------------------------*/
void weitap(sf_complex **tt  /* [n2][n1] tapered array (in and out) */,
	    weitap3d tap)
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

