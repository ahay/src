/* Non-stationary prediction-error filter estimation */
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

#include "npef.h"
#include "nhconest.h"
#include "npolydiv2.h"

#include "nhelix.h"
/*^*/

void nfind_pef(int nd     /* data size */, 
	       float *dd  /* data */, 
	       nfilter aa /* estimated filter */, 
	       nfilter rr /* regularization filter */, 
	       int niter  /* number of iterations */, 
	       float eps  /* regularization parameter */, 
	       int nh     /* filter size */) 
/*< estimate non-stationary PEF >*/
{
    int ip, ih, na, np, nr;
    float *flt;

    np = aa->np;
    nr = np*nh;
    flt = sf_floatalloc(nr);

    nhconest_init(dd, aa, nh);
    npolydiv2_init( nr, rr);

    sf_solver_prec(nhconest_lop, sf_cgstep, npolydiv2_lop,
		   nr, nr, nd, flt, dd, niter, eps, "end");
    sf_cgstep_close();
    npolydiv2_close();

    for (ip=0; ip < np; ip++) {
	na = aa->hlx[ip]->nh;
	for (ih=0; ih < na; ih++) {
	    aa->hlx[ip]->flt[ih] = -flt[ip*nh + ih];
	}
    }

    free(flt);
}
