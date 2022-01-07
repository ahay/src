/* Complex Causal integration */
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

#include "ccausint.h"

void sf_ccausint_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i;       
    sf_complex t;

    sf_cadjnull (adj, add, nx, ny, xx, yy);

    t = sf_cmplx(0.,0.);
    if (adj) {
	for (i=nx-1; i >= 0; i--) {
    #ifdef SF_HAS_COMPLEX_H
	    t += yy[i];
	    xx[i] += t;
    #else
      t = sf_cadd(t,yy[i]);
      xx[i] = sf_cadd(xx[i],t);
    #endif

	}
    } else {
	for (i=0; i <= nx-1; i++) {
    #ifdef SF_HAS_COMPLEX_H
	    t += xx[i];
	    yy[i] += t;
    #else
      t = sf_cadd(t, xx[i]);
      yy[i] = sf_cadd(yy[i],t);
    #endif

	}
    }
}

/* 	$Id$	 */
