/* Solve data space recovery using nonstationary matching radon operator. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "nmrt.h"
#include "matching.h"
#include "radonoper.h"
#include "slant.h"

static float *tmp1, *tmp2;
static int nt, nx, np;
static bool freq;

void nmrt_init (int nt1, int nx1, 
		int nf1,
		int np1,
		float dt1, float t01,
		float dx1, float ox1, float x01,
		float dp1, float p01,
		bool par1,
		bool freq1,
		bool rho1,
		float anti1,
		float p11,
		float *filt)
/*< initialize >*/
{
    nt = nt1;
    nx = nx1;
    np = np1;
    freq = freq1;
    
    tmp1 = sf_floatalloc(nt*nx);
    tmp2 = sf_floatalloc(nt*np);
   
    matching_init(filt, nt1, nx1, nf1);

    if (freq) {
	radonoper_init (nt1,dt1,t01,nx1,dx1,ox1,x01,np1,dp1,p01,par1);
    } else {
	slant_init (true,rho1,x01,dx1,nx1,p01,dp1,np1,t01,dt1,nt1,p11,anti1);
    }
}

void nmrt_lop (int niter         /* number of iterations */, 
	       float *xx        /* model */,
	       float *yy         /* data */, 
	       bool verb         /* verbosity flag */) 
/*< least squares solver >*/
{
    sf_solver (drecov_oper, sf_cgstep, nt*nx, nt*nx, xx, yy, niter, 
	       "x0", yy, "verb", verb, "end");
    sf_cgstep_close();
}

void drecov_oper (bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< linear operator >*/
{
    sf_adjnull (adj, add, nm, nd, mod, dat);

    if (adj) {
	if (freq) {
	    radonoper_lop (true, false, nt*np, nt*nx, tmp2, dat);
	    radonoper_lop (false, add, nt*np, nt*nx, tmp2, tmp1);
	} else {
	    slant_lop(true,false,nt*np,nt*nx,tmp2,dat);	
	    slant_lop(false,add,nt*np,nt*nx,tmp2,tmp1);
	}
	pmatching_lop(true,false,nt*nx,nt*nx,mod,tmp1);
    } else {
	pmatching_lop(false,false,nt*nx,nt*nx,mod,tmp1);
	if (freq) {
	    radonoper_lop (true, false, nt*np, nt*nx, tmp2, tmp1);
	    radonoper_lop (false, add, nt*np, nt*nx, tmp2, dat);
	} else {
	    slant_lop(true,false,nt*np,nt*nx,tmp2,tmp1);
	    slant_lop(false,add,nt*np,nt*nx,tmp2,dat);
	}
    }
}

void nmrt_oper (bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< linear operator >*/
{
    sf_adjnull ((bool)!adj, add, nm, nd, mod, dat);

    if (adj) {
	if (freq) {
	    radonoper_lop (false, add, nt*np, nt*nx, mod, tmp1);
	} else {
	    slant_lop(false,add,nt*np,nt*nx,mod,tmp1);
	}
	pmatching_lop(true,false,nt*nx,nt*nx,dat,tmp1);
    } else {
	pmatching_lop(false,false,nt*nx,nt*nx,dat,tmp1);

	if (freq) {
	    radonoper_lop (true, add, nt*np, nt*nx, mod, tmp1);
	} else {
	    slant_lop(true, add, nt*np, nt*nx, mod, tmp1);
	}
    }

}

void nmrt_close (void)
/*< free allocated storage >*/
{
    free (tmp1);
    free (tmp2);
}
/* 	$Id: nmrt.c 5890 2010-05-03 21:13:12Z yang_liu $	 */

