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
#include "helix_tcaf.h"

#include "match.h"
#include "copy_scaled.h"

void find_match(int nd       /* data size */,
	      float* mm        /* match [nd] */,
		  float* dd        /* data [nd] */,
	      sf_filter aa     /* estimated filter */,
	      int niter        /* number of iterations */,
	      float eps        /* regularization parameter */)
/*< find matched filter >*/
{
	int id, nh2;
	float *dd_pad=NULL,*tmp=NULL;
	float eps_S=0.0;

	bool CG = true;

	nh2 = (aa->nh-1) / 2;
	dd_pad = sf_floatalloc(nd+aa->nh-1);

	for (id=0;id<(nd+aa->nh-1);id++) {
		if (CG==true)
			dd_pad[id] = 0.0;
		else {
			dd_pad[id] = eps;
			if (id<nd) mm[id]+=eps;
		}
	}

	for (id=0;id<nd;id++)
		dd_pad[id+nh2] += dd[id];

    if (CG) {
    /* solution using preconditioned CG */ //FIXME: still memory leaks in sf_cdstep (sf_list) not affordable with huge datasets
    helix_tcaf_init(nd, mm, aa);
    sf_cdstep_init();

    sf_solver_prec(helix_tcaf_lop, sf_cdstep, sf_copy_lop, aa->nh,aa->nh, nd+aa->nh-1, aa->flt, dd_pad,
    	      niter*3, eps, "x0", aa->flt, "verb",false,"end");

    sf_cdstep_close();
    }
    else {
    /* solution using shaping CG with S = (1+eps^2)^(-1) I*/

    eps_S = 1.0 / (1.0+eps);

    tmp = sf_floatalloc(aa->nh);

    helix_tcaf_init(nd,mm, aa);
    copy_scaled_init(eps_S);

    sf_conjgrad_init(aa->nh, aa->nh, nd+aa->nh-1, nd+aa->nh-1, 1.0, 1.e-8, false, false);

	sf_conjgrad(NULL,helix_tcaf_lop,copy_scaled_lop,
			tmp,aa->flt,dd_pad,3*niter);

    sf_conjgrad_close();
    free(tmp);
    free(dd_pad);
    }

}

/* 	$Id: pef.c 2521 2007-02-02 00:25:42Z sfomel $	 */

