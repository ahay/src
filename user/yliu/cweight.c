/* Weight operator for complex numbers. */
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
/*^*/

#include "cweight.h"

static float* w;

void cweight_init(float *w1)
/*< initialize weight >*/
{
    w = w1;
}

void cweight_lop (bool adj, bool add, int nx, int ny, 
		 sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	if (adj) {
	    xx[i] += yy[i] * w[i];
	} else {
	    yy[i] += xx[i] * w[i];
	}
#else
	if (adj) {
	    xx[i] = sf_cadd(xx[i],sf_crmul(yy[i],w[i]));
	} else {
	    yy[i] = sf_cadd(yy[i],sf_crmul(xx[i],w[i]));
	}
#endif
    }
}

/* 	$Id: cweight.c 1836 2006-05-09 10:38:59Z fomels $	 */
