/* Plane prediction filter on a helix. */
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

#include "wilson.h"
#include "compress.h"

static int nt, nx;
static sf_filter ss;

void pfactor_init(int nt1, int nx1 /* data size */)
/*< initialize >*/
{
    nt = nt1;
    nx = nx1;

    ss =  sf_allocatehelix(8);
    ss->lag[0] = 1;
    ss->lag[1] = 2;
    ss->lag[2] = nt - 1;
    ss->lag[3] = nt;
    ss->lag[4] = nt + 1;
    ss->lag[5] = nt * nx - 1;
    ss->lag[6] = nt * nx;
    ss->lag[7] = nt * nx + 1;

    wilson_init (nt*nx*10);
}

void pfactor_close(void)
/*< free allocated storage >*/
{
    free(ss);
    wilson_close();
}

sf_filter pfactor(int na     /* filter size */,
		  float p    /* inline slope */,
		  float q    /* crossline slope */,
		  int niter  /* number of iterations */,
		  float eps  /* compression tolerance */,
		  bool fixed /* if fixed size */,
		  bool verb  /* verbosity */,
		  float *a0  /* zero-lag coefficient */)
/*< find a filter >*/
{
    sf_filter aa;
    int i;
    float f, s0;

    f = p*p*(p*p-1.) + q*q*(q*q-1.);
    s0 = 2. + 0.75*f;
    ss->flt[0] = -f;
    ss->flt[1] = 0.25*f;
    ss->flt[2] = 0.5*p*(1.-p);
    ss->flt[3] = p*p-1.;
    ss->flt[4] = -0.5*p*(1.+p);
    ss->flt[5] = 0.5*q*(1.-q);
    ss->flt[6] = q*q-1.;
    ss->flt[7] = -0.5*q*(1.+q);

    if (fixed) {
	na = 25;
	aa = sf_allocatehelix(na);

	aa->lag[0] = 1;
	aa->lag[1] = 2;
	aa->lag[2] = nt-1;
	aa->lag[3] = nt;
	aa->lag[4] = nt+1;
	aa->lag[5] = 2*nt-1;
	aa->lag[6] = 2*nt;
	aa->lag[7] = 2*nt+1;
	aa->lag[8]  = nt*(nx-3)-2;
	aa->lag[9] = nt*(nx-3)-1;
	aa->lag[10] = nt*(nx-3);
	aa->lag[11] = nt*(nx-3)+1;
	aa->lag[12] = nt*(nx-3)+2;
	aa->lag[13] = nt*(nx-2)-2;
	aa->lag[14] = nt*(nx-2)-1;
	aa->lag[15] = nt*(nx-2);
	aa->lag[16] = nt*(nx-2)+1;
	aa->lag[17] = nt*(nx-2)+2;
	aa->lag[18] = nt*(nx-1)-1;
	aa->lag[19] = nt*(nx-1);
	aa->lag[20] = nt*(nx-1)+1;
	aa->lag[21] = nt*(nx-1)+2;
	aa->lag[22] = nt*nx-1;
	aa->lag[23] = nt*nx;
	aa->lag[24] = nt*nx+1;
    } else {
	aa = sf_allocatehelix(na);
	for (i=0; i < na; i++) {
	    aa->lag[i] = i+1;
	}
    }

    for (i=0; i < na; i++) {
	aa->flt[i] = 0.;
    }

    *a0 =  wilson_factor (niter, 2.*s0, ss, aa, verb, 1.e-6);
    aa = compress(aa,eps);

    if (!fixed) {
	na = aa->nh;
	for (i=0; i < na; i++) {
	    aa->flt[i] = 0.;
	}

	*a0 =  wilson_factor (niter, 2.*s0, ss, aa, verb, 1.e-6);
    }

    return aa;
}
