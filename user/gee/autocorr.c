/* Auto-correlation of a helix filter. */
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

#include "autocorr.h"
#include "compress.h"

sf_filter autocorr(const sf_filter aa /* input filter */, 
		   float a0        /* input zero lag */, 
		   float *s0       /* output zero lag */, 
		   float eps       /* tolerance for compression */)
/*< Output the autocorrelation (positive side only) of the input filter >*/
{
    int i, j, k, n, na;
    float f, b0;
    sf_filter ss;

    na = aa->nh;

    ss = sf_allocatehelix (na*(na+1)/2);
 
    b0 = a0*a0;
    for (i=0; i < na; i++) {
	f = aa->flt[i];
	b0 += f*f;
	ss->flt[i] = a0*f;
	ss->lag[i] = aa->lag[i];
    }
    *s0 = b0;
	  
    k = na-1;
    for (i=0; i < na; i++) {
	for (j = i+1; j < na; j++) {
	    k++;

	    ss->flt[k] = aa->flt[j] * aa->flt[i];
	    ss->lag[k] = aa->lag[j] - aa->lag[i];

	    for (n=0; n < k-1; n++) {
		if (ss->lag[n] == ss->lag[k] && ss->flt[n] != 0.) {
		    ss->flt[n] += ss->flt[k];
		    ss->flt[k] = 0.;
		}
	    }
	}
    }

    return compress(ss,eps);
}
