/* RTM of incomplete zero-offset data with LS and shaping regularization */
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

#include "postrtm.h"
#include "pwdsl.h"

static float *tmp;

void ptshape_init(int nx, int nz, int nt, bool verb, int n0, int padx, int padz, 
                 int padnx, int padnz, float dx, float dz, int *head, float **padvv, 
                    int rect1, int rect2, int order, float **dp, float lambda)
/*< initialize >*/
{
    int nxz, nxt;
    
    nxz=nx*nz;
    nxt=nx*nt;
    
    postrtm_init(nx, nz, nt, n0, padx, padz, padnx, padnz, dx, dz, head, padvv);
    pwdsl_init(nz, nx, order, rect1, rect2, 0.01);
    pwdsl_set(dp);
    
    tmp=sf_floatalloc(nxz);
    sf_conjgrad_init(nxz, nxz, nxt, nxt, lambda, 10*FLT_EPSILON, verb, false);
}

void ptshape_close()
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    pwdsl_close();
    postrtm_close();
}

void ptshape(int niter, int nxz, int nxt, float *mm, float *dd)
/*< conjugate-gradient with shaping >*/
{
    sf_conjgrad(NULL, postrtm_lop, pwdsl_lop, tmp, mm, dd, niter);
}
