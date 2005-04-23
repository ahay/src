/* 2-D inverse interpolation with PEF estimation */
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

#include "levint2.h"

#include "helix.h"
/*^*/

#include "helicon.h"
#include "pefhel.h"
#include "lint2.h"

void levint2 (bool cdstep                /* if use conj. directions */,
	      int niter                  /* number of iterations */, 
	      int warmup                 /* initial iterations */,
	      int nd                     /* data size */,
	      float *x, float *y         /* data coordinates */, 
	      float *dd                  /* data */, 
	      float o1, float d1, int n1 /* inline */, 
	      float o2, float d2, int n2 /* crossline */, 
	      filter aa                  /* PEF */, 
	      float *rr                  /* residual */, 
	      float eps                  /* regularization parameter */) 
/*< 2-D inverse interpolation with PEF estimation >*/
{
    int nm, nr;

    nm = n1*n2;
    nr = nm+aa->nh;

    lint2_init (n1,o1,d1, n2,o2,d2, x, y);
    pefhel_init (aa, nm, rr);

    if (cdstep) {
	sf_cdstep_init();
	sf_solver_reg (lint2_lop, sf_cdstep, helicon_lop, nm, nm, nd, rr, dd, 
		       warmup, eps, "x0", rr, "nmem", warmup, "end");
	sf_cdstep_close();
	/*
	  sf_solver_reg (lint2_lop, sf_cdstep, pefhel_lop, nm, nr, nd, rr, dd, 
	  niter, eps, "x0", rr, "nlreg", helicon_lop, "nmem", 3, "end");
	  sf_cdstep_close();
	*/
    } else {
	sf_solver_reg (lint2_lop, sf_cgstep, helicon_lop, nm, nm, nd, rr, dd, 
		       warmup, eps, "x0", rr, "end");
	sf_cgstep_close();
	/*
	  sf_solver_reg (lint2_lop, sf_cgstep, pefhel_lop, nm, nr, nd, rr, dd, 
	  niter, eps, "x0", rr, "nlreg", helicon_lop,"end");
	  sf_cgstep_close();
	*/
    }
}
