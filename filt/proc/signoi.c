/* Signal and noise separation with helical prediction-error filters */
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
/*^*/

#include "signoi.h"

#include "helicon.h"
#include "polydiv.h"

#include "helix.h"
/*^*/

static int niter, nd;
static float eps, *dd;
static filter nn, ss;


void signoi_init(filter nn_in /* noise PEF */, 
		 filter ss_in /* signal PEF */, 
		 int niter_in /* number of iterations */, 
		 int nd_in    /* data size */, 
		 float eps_in /* regularization parameter (signal/noise) */)
/*< initialize >*/
{
    nn = nn_in;
    ss = ss_in;
    niter = niter_in;
    nd = nd_in;
    eps = eps_in;
    dd = sf_floatalloc(nd);
}

void signoi_lop (bool adj, bool add, int n1, int n2, float *data, float *sign)
/*< linear operator >*/
{
    helicon_init (nn);
    polydiv_init (nd, ss); 

    if (nd != n1 || nd != n2) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,n1,n2,data,sign);

    helicon_lop (false, false, n1, n1, data, dd);
    sf_solver_prec(helicon_lop, sf_cgstep, polydiv_lop, nd, nd, nd, sign, dd,
		   niter, eps, "end");
    sf_cgstep_close();

    nn++;
    ss++;
}

void signoi_close(void)
/*< free allocated storage >*/
{
    free(dd);
    polydiv_close();
}
