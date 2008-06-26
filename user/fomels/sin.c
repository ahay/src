/* Operations with complex sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "sin.h"
#include "ccopy.h"

static sf_complex z0;

void sin_init(sf_complex z1)
/*< initialize >*/
{
    z0 = z1;
}

void sin_destruct(bool adj, bool add, int nx, int ny, 
		  sf_complex *xx, sf_complex *yy)
/*< destruction operator >*/
{
    int i;

    sf_cadjnull(adj, add, nx, nx, xx, yy);
    
    for (i = 1; i < nx; i++) {
	if(adj) {
#ifdef SF_HAS_COMPLEX_H	 
	    xx[i]   += yy[i];
	    xx[i-1] -= yy[i] * conjf(z0);
#else
	    xx[i-1] = sf_cadd(xx[i-1],sf_cmul(yy[i],sf_cneg(sf_conjf(z0))));
#endif
	} else {
#ifdef SF_HAS_COMPLEX_H	 
	    yy[i] += xx[i] - xx[i-1] * z0;
#else
	    yy[i] = sf_cadd(yy[i],sf_cmul(xx[i-1],sf_cneg(z0)));
#endif
	}
    }
}
