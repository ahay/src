/* Matched filter estimation */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include "helix_icaf.h"

#include "match.h"
#include "copy_scaled.h"

void find_match(int nd       /* data size */,
	      float* mm        /* match [nd] */,
		  float* dd        /* data [nd] */,
	      sf_filter aa     /* estimated filter */,
	      int liter        /* number of CG iterations */,
	      float eps        /* regularization parameter */)
/*< find matched filter >*/
{
    bool CG = false;

    if (CG) {
    /* solution using preconditioned CG */ //FIXME: still memory leaks in sf_cdstep (sf_list) not affordable with huge datasets
    helix_icaf_init(mm, aa, 1);
    sf_cdstep_init();

    sf_solver_prec(helix_icaf_lop, sf_cdstep, sf_copy_lop, aa->nh,aa->nh, nd, aa->flt, dd,
    	      liter*3, eps, "x0", aa->flt, "verb",false,"end");

    sf_cdstep_close();
//    sf_solver_prec(helix_icaf_lop, sf_cgstep, sf_copy_lop, aa->nh,aa->nh, nd, aa->flt, dd,
//    		liter*3, eps, "x0", aa->flt, "verb",false,"end");
    }
    else {
    /* solution using shaping CG with S = (1+eps^2)^(-1) I*/

    float *tmp;
    float eps_S = 1.0 / (1.0+eps);
    int i;

    tmp = sf_floatalloc(aa->nh);

    for (i=0;i<nd;i++) { // this is a slight dump to avoid 0 = aa->flt * 0
    	mm[i]+=eps;
    	dd[i]+=eps;
    }

    helix_icaf_init(mm, aa, 1);
    copy_scaled_init(eps_S);

    sf_conjgrad_init(aa->nh, aa->nh, nd, nd, 1.0, 1.e-8, false, false);

	sf_conjgrad(NULL,helix_icaf_lop,copy_scaled_lop,
			tmp,aa->flt,dd,3*liter);

    sf_conjgrad_close();
    free(tmp);
    }

}


