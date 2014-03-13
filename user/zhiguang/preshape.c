/* 2D prestack RTM with LS and shaping regularization */
/*
 Copyright (C) 2014 University of Texas at Austin
 
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

#include "prertm.h"
#include "pwdsl.h"

static float *tmp;

/* preshape_init(verb, nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay,
                  dx, dz, padx, padz, padnx, padnz,
                  nm, nd, dr_v, ds_v, r0_v, s0_v, zr_v, zs_v,
                  padvv, ww, rect1, rect2, order, dp, lambda); */
                  
void preshape_init(bool verb, int nz, int nx, int nt, int nr, int ns, int nw,
                   int nsource, int dsource, int ndelay, float dx, float dz,
                   int padx, int padz, int padnx, int padnz, int nm, int nd,
                   int dr_v, int ds_v, int r0_v, int s0_v, int zr_v, int zs_v,
                   float **padvv, float *ww, int rect1,
                   int rect2, int order, float **dp, float lambda)
/*< initialize >*/
{
    prertm_init(verb, nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay, dx, dz,
    padx, padz, padnx, padnz, dr_v, ds_v, r0_v, s0_v, zr_v, zs_v, padvv, ww);
    pwdsl_init(nz, nx, order, rect1, rect2, 0.01);
    pwdsl_set(dp);
    
    tmp=sf_floatalloc(nm);
    sf_conjgrad_init(nm, nm, nd, nd, lambda, 10*FLT_EPSILON, verb, false);
}

void preshape_close()
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    pwdsl_close();
    prertm_close();
}

void preshape(int niter, float *mm, float *dd)
/*< conjugate-gradient with shaping >*/
{
    sf_conjgrad(NULL, prertm2_oper, pwdsl_lop, tmp, mm, dd, niter);
}
