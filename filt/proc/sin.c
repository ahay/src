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
    z0 = conjf(z1);
}

void sin_destruct(bool adj, bool add, int nx, int ny, 
		  sf_complex *xx, sf_complex *yy)
/*< destruction operator >*/
{
    int i;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    ccopy_lop(adj,add,nx,nx,xx,yy);

    for (i = 1; i < nx; i++) {
	if(adj) {
#ifdef SF_HAS_COMPLEX_H	 
	    xx[i-1] -= yy[i] * conjf(z0);
#else
	    xx[i-1] = sf_cadd(xx[i-1],sf_cmul(yy[i],sf_cneg(sf_conjf(z0))));
#endif
	} else {
#ifdef SF_HAS_COMPLEX_H	 
	    yy[i] -= xx[i-1] * z0;
#else
	    yy[i] = sf_cadd(yy[i],sf_cmul(xx[i-1],sf_cneg(z0)));
#endif
	}
    }
}

void sin_construct(bool adj, bool add, int nx, int ny, 
		   sf_complex *xx, sf_complex *yy)
/*< construction operator >*/
{
    int i;
    sf_complex t;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    sf_cadjnull(adj, add, nx, nx, xx, yy);

    t = sf_cmplx(0.0,0.0);
    if(adj) {	
	for (i = nx-1; i >= 0; i--) {
#ifdef SF_HAS_COMPLEX_H
	    t = yy[i] + t * conjf(z0);
	    xx[i] += t;  
#else
	    t = sf_cadd(yy[i],sf_cmul(t,sf_conjf(z0)));
	    xx[i] = sf_cadd(xx[i],t);
#endif
	}
    } else {
	for (i = 0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H	    
	    t = xx[i] + t * z0;
	    yy[i] += t;
#else
	    t = sf_cadd(xx[i],sf_cmul(t,z0));
	    yy[i] = sf_cadd(yy[i],t);
#endif
	}
    }
}
