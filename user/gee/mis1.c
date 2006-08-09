/* 1-D missing data interpolation with known filter */
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

#include "mis1.h"

#include "tcai1.h"

static int nx, ny;
static float *zero;

void mis1_init(int n1    /* data length */, 
	       int na    /* filter length */, 
	       float *aa /* filter [na] */)
/*< initialize >*/
{
    int iy;

    nx=n1;
    ny=nx+na-1;

    zero = sf_floatalloc(ny);
    for (iy=0; iy < ny; iy++) {
	zero[iy]=0.;
    }

    tcai1_init(na,aa);
}

void mis1_close(void)
/*< free allocated storage >*/
{
    free(zero);
}

void mis1(int niter         /* number of iterations */, 
	  float *xx         /* data/model */, 
	  const bool *known /* mask for known data */,
	  const char *step  /* solver */) 
/*< interpolate >*/
{
    switch (step[1]) {
	case 'g': /* conjugate gradients */
	    sf_solver (tcai1_lop, sf_cgstep, nx, ny, xx, zero, 
		       niter, "x0", xx, "known", known, "end");
	    sf_cgstep_close();
	    break;
	case 'd': /* conjugate directions */
	    sf_cdstep_init();
	    sf_solver (tcai1_lop, sf_cdstep, nx, ny, xx, zero, 
		       niter, "x0", xx, "known", known, "end");
	    sf_cdstep_close();
	    break;
	default:
	    sf_error("%s: unknown step %s",__FILE__,step);
	    break;
    }
}

/* 	$Id$	 */
