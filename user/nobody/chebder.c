/* Taking derivative in the Chebyshev transform domain. */
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

void chebder (bool inv       /* inverse flag (derivative or integral) */, 
	      int n          /* data size */, 
	      const float *f /* input data */, 
	      float *b       /* output derivative */)
/*< Compute derivative or integral in the Chebyshev domain >*/
{
    int i;

    n--;
    if (inv) {
	b[n] = f[n-1]/(2.*n);
	for (i=2; i < n; i++) {
	    b[i] = (f[i-1] - f[i+1])/(2.*(i-1));
	}
	b[1] = f[0] - 0.5*f[2];
	b[0] = 0.;
	for (i=1; i < n; i+= 2) {
	    b[0] += b[i] - b[i+1];
	}
    } else {
	b[n] = 0.;
	b[n-1] = 2.*n*f[n];
	for (i=n-1; i >= 2; i--) {
	    b[i-1] = 2.*(i-1)*f[i] + b[i+1];
	}
	b[0] = f[1] + 0.5*b[2];
    }
}

