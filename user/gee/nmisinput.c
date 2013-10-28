/* Missing data mask for a non-stationary helix filter */
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

#include "nmisinput.h"
#include "nhelicon.h"

#include "nhelix.h"
/*^*/

void nfind_mask(int nd           /* data size */, 
		const int *known /* mask for known data [nd] */, 
		nfilter aa       /* filter */) 
/*< find mask >*/
{
    float *rr, *dfre;
    int ip, i;

    rr = sf_floatalloc(nd);
    dfre = sf_floatalloc(nd);

    for (i=0; i < nd; i++) {
	dfre[i] = known[i]? 0.:1.;
    }

    nhelicon_init( aa);

    for (ip=0; ip < aa->np; ip++) {
	for (i=0; i < aa->hlx[ip]->nh; i++) {
	    aa->hlx[ip]->flt[i] = 1.;
	}
    }
	
    
    nhelicon_lop(false,false,nd,nd,dfre,rr);
    for (i=0; i < nd; i++) {
	if ( rr[i] > 0.) aa->mis[i] = true;
    }

    for (ip=0; ip < aa->np; ip++) {
	for (i=0; i < aa->hlx[ip]->nh; i++) {
	    aa->hlx[ip]->flt[i] = 0.;
	}
    }
}
