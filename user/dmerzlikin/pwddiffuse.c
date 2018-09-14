/* Anisotropic diffusion by regularized inversion. PWD instead of gradients. */
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
#include <assert.h>

#include <rsf.h>
/*^*/
#include "allp3.h"

static float *vx, *vy, eps, *pwdx, *pwdy, *pwddata, *pp1, *pp2;
static int niter, repeat;
static int n1, n2, n3, n123;
static allpass ap, aq;

void pwdgradient (bool adj, bool add, 
		      int nx, int ng, float* x, float* g)
/*< gradient operator not for external use >*/
{
    int i;

    assert(nx == n123 && n123 == n3*n2*n1 && ng == n123);

    sf_adjnull (adj,add,nx,ng,x,g);

    if (adj == false){/* forward */

	allpass32_lop(adj,false,n123,n123,x,pwdx,x,pwdy);

    	for (i=0; i < n123; i++){

		g[i] += vx[i]*pwdx[i] + vy[i]*pwdy[i];

	}

    } else {/* adjoint */
	
	for (i=0; i < n123; i++){

		pwddata[i] = g[i];
		g[i] = vx[i]*g[i];
		pwddata[i] = vy[i]*pwddata[i];

	}

	allpass32_lop(adj,false,n123,n123,pwdx,g,pwdy,pwddata);


	for (i=0; i < n123; i++){

		x[i] += pwdx[i] + pwdy[i];

	}

    }

}

void pwddiffuse_init(int nt, int nx, int ny /* data size */,
		 float *fpp1, float *fpp2 /* inline and crossline dips */, 
		 int nw, int nj1, int nj2 /* PWD parameters */,
		 float *fvx, float *fvy /* parallel to edges vector components */,
		 int fniter /* number of iterations */,
		 int frepeat /* number of smoothing iterations */,
		 float feps /* regularization */)
/*< initialize anisotropic diffusion >*/
{

    vx = fvx;
    vy = fvy;

    pp1 = fpp1;
    pp2 = fpp2;

    niter = fniter;
    repeat = frepeat;
  
    eps = feps;

    n1 = nt;
    n2 = nx;
    n3 = ny;
    n123 = n1*n2*n3;

    /* iline */
    ap = allpass_init(nw, nj1, n1,n2,n3, pp1);

    /* xline */
    aq = allpass_init(nw, nj2, n1,n2,n3, pp2);

    /* iline + xline */
    allpass3_init(ap,aq);

    pwdx = sf_floatalloc(n123);
    pwdy = sf_floatalloc(n123);
    pwddata = sf_floatalloc(n123);

}

void pwddiffuse_fwdoradj(bool adj, bool add, int fnx, int fny, float* x, float* y)
/*< anisotropic diffusion loop >*/
{

    int i;

    if((n123 != fnx) || (n123 != fny)) sf_error("Wrong dimensions");

    sf_adjnull(adj,add,fnx,fny,x,y);

    for (i=0; i < n123; i++){

	pwdx[i] = 0.0;
	pwdy[i] = 0.0;
	pwddata[i] = 0.0;

    }

    pwdgradient(adj,add,fnx,fny,x,y);

}

void pwddiffuse_lop(int fnx, int fny, float* x, float* y)
/*< anisotropic diffusion loop >*/
{

    int i;

    if((n123 != fnx) || (n123 != fny)) sf_error("Wrong dimensions");

    for (i=0; i < n123; i++){

	pwdx[i] = 0.0;
	pwdy[i] = 0.0;
	pwddata[i] = 0.0;

    }

    /*for (i1 = 0; i1 < n123; i1++){

	y[i1] = 0.0;

    }*/


    //sf_warning("before repeat");

    for (i=0; i < repeat; i++) {
	/* Call a generic function for solving
	   .   F[m] =~ d
	   eps*D[m] =~ 0 
	*/

	sf_solver_reg(sf_copy_lop /* F */,
		      sf_cgstep   /* CG step */,
		      pwdgradient    /* D */,
		      n123         /* size of D[m] */,
		      n123         /* size of m */,
		      n123         /* size of d */,
		      y        /* m (output) */,
		      x        /* d (input)  */,
		      niter       /* # iterations */,
		      eps         /* eps */,
		      "verb",false,"end");

	sf_cgstep_close(); /* free internal storage */

    }/* repeat loop */

}

void pwddiffuse_close(void)
/*< free allocated storage >*/
{
    free(pwdx);
    free(pwdy);
    free(pwddata);
    sf_warning("pwddiffuse: memory free");
}
