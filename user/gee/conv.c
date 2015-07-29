/* Convolution of helix filters */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "compress.h"

sf_filter conv(const sf_filter aa, 
	       const sf_filter bb,
               bool one /* include leading one */) 
/*< convolve aa and bb >*/
{
    sf_filter ss;
    int na, nb, amax, bmax, ns, ia, ib, i;
        
    na = aa->nh;
    nb = bb->nh;

    amax = 0;
    for (i=0; i < na; i++) {
	if (aa->lag[i] > amax) amax = aa->lag[i];
    }
    bmax = 0;
    for (i=0; i < nb; i++) {
	if (bb->lag[i] > bmax) bmax = bb->lag[i];
    }
    ns = amax+bmax;
    ss = sf_allocatehelix(ns);
    for (i=0; i < ns; i++) {
	ss->lag[i] = i+1;
	ss->flt[i] = 0.;
    }

    if (one) {
	for (ia=0; ia < na; ia++) {
	    i = aa->lag[ia];
	    ss->flt[i-1] = aa->flt[ia];
	}
	for (ib=0; ib < nb; ib++) {
	    i = bb->lag[ib];
	    ss->flt[i-1] += bb->flt[ib];
	}
    }

    for (ia=0; ia < na; ia++) {
	for (ib=0; ib < nb; ib++) {
	    i = aa->lag[ia] + bb->lag[ib]; 
	    ss->flt[i-1] += aa->flt[ia] * bb->flt[ib];
	}
    }

    return compress(ss,FLT_EPSILON);
}

