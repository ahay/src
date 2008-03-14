/* Fill in missing data using non-stationary filter. */
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

#include "nhelix.h"
/*^*/

#include "nhelicon.h"
#include "npolydiv.h"

void nmis2( int niter        /* number of iterations */, 
	    int nx           /* data size */,
	    float *xx        /* in - data with holes, out - filled */, 
	    const nfilter aa /* filter */, 
	    bool *known      /* mask for known data */, 
	    bool precon      /* whether to use preconditioning */) 
/*< fill missing data >*/
{
    int ix;
    float *dd;

    if (precon) {  /*  KP p = K y,   m = P p */
	sf_mask_init( known);
	npolydiv_init( nx, aa);

	sf_solver_prec( sf_mask_lop, sf_cgstep, npolydiv_lop,
			nx, nx, nx, xx, xx, niter, 0., "end");
	npolydiv_close();
    } else {       /*  KA m = 0 */
	dd = sf_floatalloc(nx);
	for (ix=0; ix < nx; ix++) {
	    dd[ix] = 0.;
	}

	nhelicon_init( aa);
	sf_solver( nhelicon_lop, sf_cgstep, 
		   nx, nx, xx, dd, niter, 
		   "known", known, "x0", xx, "end");
	free(dd);
    }
    sf_cgstep_close();
}

