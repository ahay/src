/* Simple weight operator */
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

#include "_bool.h"
#include "error.h"
#include "adjnull.h"
#include "komplex.h"

#include "weight.h"

static float* w;

void sf_weight_init(float *w1)
/*< initialize >*/
{
    w = w1;
}

void sf_weight_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
        for (i=0; i < nx; i++) {
	    xx[i] += yy[i] * w[i];
	}
    } else {
        for (i=0; i < nx; i++) {
            yy[i] += xx[i] * w[i];
	}
    }

}

void sf_cweight_lop (bool adj, bool add, int nx, int ny, 
		     sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
        for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    xx[i] += yy[i] * w[i];
#else
	    xx[i] = sf_cadd(xx[i],sf_crmul(yy[i],w[i]));
#endif
        }
    } else {
        for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    yy[i] += xx[i] * w[i];
#else
	    yy[i] = sf_cadd(yy[i],sf_crmul(xx[i],w[i]));
#endif
        }
    }

}

void sf_weight_apply(int nx, float *xx)
/*< apply weighting in place >*/
{
    int i;

    for (i=0; i < nx; i++) {
	xx[i] *= w[i]*w[i];
    }
}


void sf_cweight_apply(int nx, sf_complex *xx)
/*< apply weighting in place >*/
{
    int i;

    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
      xx[i] *= w[i]*w[i];
#else
      xx[i] = sf_crmul(xx[i],w[i]*w[i]);
#endif
    }
}

/* 	$Id$	 */
