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
#include "helicon2.h"

static bool verb;
static int niter, nd;
static float eps, *dd;
static sf_filter nn, ss;


void signoi_init(sf_filter nn_in /* noise PEF */, 
		 sf_filter ss_in /* signal PEF */, 
		 int niter_in    /* number of iterations */, 
		 int nd_in       /* data size */, 
		 float eps_in    /* regularization parameter (signal/noise) */,
		 bool verb_in    /* verbosity flag */)
/*< initialize >*/
{
    nn = nn_in;
    ss = ss_in;
    niter = niter_in;
    nd = nd_in;
    eps = eps_in;
    verb = verb_in;
    dd = sf_floatalloc(nd);
}

void signoi_lop (bool adj, bool add, int n1, int n2, 
		 float *data, float *sign)
/*< linear operator >*/
{
    sf_helicon_init (nn);
    sf_polydiv_init (nd, ss); 

    sf_adjnull(adj,add,n1,n2,data,sign);

    sf_helicon_lop (false, false, n1, n1, data, dd);
    sf_solver_prec(sf_helicon_lop, sf_cgstep, sf_polydiv_lop, 
		   nd, nd, nd, sign, dd, niter, eps, 
		   "verb", verb, "end");
    sf_cgstep_close();

    nn++;
    ss++;
}

void signoi2_lop (bool adj, bool add, int n1, int n2, 
		 float *data, float *sign)
/*< alternative linear operator >*/
{
    helicon2_init (nd, nn);
    sf_helicon_init (ss); 

    sf_adjnull(adj,add,n1,n2,data,sign);

    helicon2_lop (false, false, n1, n1, data, dd);
    sf_solver_reg(helicon2_lop, sf_cgstep, sf_helicon_lop, 
		  nd, nd, nd, sign, dd, niter, eps, 
		  "verb", verb, "end");
    sf_cgstep_close();

    nn++;
    ss++;
}

void signoi1_lop (bool adj, bool add, int n1, int n2, 
		 float *data, float *sign)
/*< alternative linear operator >*/
{
    helicon2_init (nd, nn);
    sf_polydiv_init (nd, ss); 

    sf_adjnull(adj,add,n1,n2,data,sign);

    helicon2_lop (false, false, n1, n1, data, dd);
    sf_solver_prec(helicon2_lop, sf_cgstep, sf_polydiv_lop, 
		   nd, nd, nd, sign, dd, niter, eps, 
		   "verb", verb, "end");
    sf_cgstep_close();
    
    nn++;
    ss++;
}

void signoi_close(void)
/*< free allocated storage >*/
{
    free(dd);
    sf_polydiv_close();
}
