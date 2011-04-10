/* Find 1-D missing data together with a prediction-error filter */
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

#include "misif.h"

#include "peftc.h"
#include "tcai2.h"

void misif1 (int niter /* number of iterations */, 
	     int na    /* filter length */, 
	     int nx    /* data length */, 
	     float *xx /* data */, 
	     float *aa /* filter */, 
	     bool *mm  /* mask for knowh data */) 
/*< interpolate >*/
{
    float *bb, *dat, *x;
    int id, nd;

    nd = nx+na;

    dat = sf_floatalloc(nd);
    x =  sf_floatalloc(nd);
    bb = x+nx;

    for (id=0; id < nd; id++) {
	dat[id] = 0.;
    }

    for (id=0; id < nx; id++) {
	x[id] = mm[id]? xx[id]: 0.;
    }
    for (id=nx; id < nd; id++) {
	x[id] = mm[id]? 1.: 0.;
    }

    peftc_init (na, nx, bb, x);
    tcai2_init (na, nx, bb);

    sf_solver (peftc_lop, sf_cgstep, nd, nd, x, dat, niter, "x0", x, 
               "nloper", tcai2_lop, "known", mm, "end");
    sf_cgstep_close ();

    for (id=0; id < nx; id++) {
	xx[id] = x[id];
    }
    for (id=nx; id < nd; id++) {
	aa[id-nx] = x[id];
    }

    free(x);
    free(dat);
}
