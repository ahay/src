/* Double integration */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "doubint.h"

#include "_bool.h"
/*^*/

void sf_doubint(bool dble, int n, float *trace /* [n] */)
/*< double integration or causal integration (in place) >*/
{
    float t;
    int it;

    
    t = 0.;
    for (it = 0; it < n; it++) {
	t += trace[it];
	trace[it] = t;
    }

    if (dble) {
	t = 0.;
	for (it = n-1; it >=0; it--) {
	    t += trace[it];
	    trace[it] = t;
	}
    }
}
