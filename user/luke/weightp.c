/* Simple weight operator with openmp */
/*
  Copyright (C) 2018 University of Texas at Austin
  
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "adjnullp.h"
/*^*/
static float* w;

void sf_weightp_init(float *w1)
/*< initialize >*/
{
    w = w1;
}

void sf_weightp_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_adjnullp (adj, add, nx, ny, xx, yy);
  
    if (adj) {
#ifdef _OPENMP
#pragma omp parallel  for
#endif
        for (i=0; i < nx; i++) {
	    xx[i] += yy[i] * w[i];
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel  for
#endif
        for (i=0; i < nx; i++) {
            yy[i] += xx[i] * w[i];
	}
    }

}

void sf_cweightp_lop (bool adj, bool add, int nx, int ny, 
		     sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    xx[i] += yy[i] * w[i];
#else
	    xx[i] = sf_cadd(xx[i],sf_crmul(yy[i],w[i]));
#endif
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    yy[i] += xx[i] * w[i];
#else
	    yy[i] = sf_cadd(yy[i],sf_crmul(xx[i],w[i]));
#endif
        }
    }

}

void sf_weightp_apply(int nx, float *xx)
/*< apply weighting in place >*/
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx; i++) {
	xx[i] *= w[i]*w[i];
    }
}


void sf_cweightp_apply(int nx, sf_complex *xx)
/*< apply weighting in place >*/
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
      xx[i] *= w[i]*w[i];
#else
      xx[i] = sf_crmul(xx[i],w[i]*w[i]);
#endif
    }
}

/* 	$Id$	 */
