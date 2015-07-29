/* Non-stationary 3-D plane-wave filter */
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
#include "pfactor.h"
#include "regrid.h"

nfilter npfactor(int npx, int npy /* number of slopes */,
		 const int *n     /* [3] data size */, 
		 const int *m     /* [3] patch size */,
		 const float *pp  /* inline slope */,
                 const float *qq  /* crossline slope */,
		 int niter        /* number of iterations */,
                 float eps        /* compression tolerance */)
/*< create a filter >*/
{
    nfilter pfilt;
    sf_filter aa;
    int i, iy, ix, ntxy, ptx;
    float pmin, pmax, qmin, qmax, dp, dq, a0, p, q;

    ntxy = n[0]*n[1]*n[2];
    ptx = m[0]*m[1];

    pfilt = (nfilter) sf_alloc(1,sizeof(*pfilt));
    pfilt->hlx = (sf_filter*) sf_alloc(npx *npy,sizeof(sf_filter));
    pfilt->pch = sf_intalloc(ntxy);
    pfilt->mis = NULL;

    pfactor_init(m[0],m[1]);

    pmin = pp[0]; pmax = pp[0];
    qmin = qq[0]; qmax = qq[0];
    for (i=1; i < ntxy; i++) {
	if (pp[i] < pmin) pmin=pp[i];
	if (qq[i] < qmin) qmin=qq[i];
	if (pp[i] > pmax) pmax=pp[i];
	if (qq[i] > qmax) qmax=qq[i];
    }
    dp = (pmax-pmin)/(npx-1);
    dq = (qmax-qmin)/(npy-1);

    for (iy=0; iy < npy; iy++) {
	q = qmin + iy * dq;
	
	for (ix=0; ix < npx; ix++) {
	    p = pmin + ix * dp;

	    aa = pfactor(ptx+1,p,q,niter,eps,false,false,&a0);
	    sf_warning("%d %d %d", ix, iy, aa->nh);

	    regrid(3, m, n, aa);

	    i = ix + npx * iy;
	    pfilt->hlx[i] = aa;
	}
    }

    pfactor_close();

    for (i=0; i < ntxy; i++) {
	ix = 0.5 + (pp[i] - pmin)/ dp;
	iy = 0.5 + (qq[i] - qmin)/ dq;
	pfilt->pch[i] = ix + npx * iy;
    }

    return pfilt;
}

