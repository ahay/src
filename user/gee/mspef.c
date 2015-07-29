/* Multi-scale prediction-error filter estimation */
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

#include "mspef.h"
#include "mshconest.h"

#include "mshelix.h"
/*^*/

void msfind_pef(int nd      /* data size */, 
		float* dd   /* data */, 
		msfilter aa /* estimated filter */, 
		int niter   /* number of iterations */) 
/*< estimate PEF >*/
{
    float *ee;
    int is, id, ns;

    ns = aa->ns;
    ee = sf_floatalloc(nd*ns);
    for (is=0; is < ns; is++) {
	for (id=0; id < nd; id++) {
	    ee[id+is*nd] = dd[id];
	}
    }

    mshconest_init( dd, aa);
    sf_solver(mshconest_lop, sf_cgstep, aa->nh, nd*ns, aa->flt, 
	      ee, niter, "x0", aa->flt, "end");
    sf_cgstep_close();

    free(ee);
}

/* 	$Id: mspef.c 11207 2013-10-28 16:56:21Z sfomel $	 */

