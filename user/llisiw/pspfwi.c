/* Full-waveform Inversion Interface */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "helm.h"
#include "pspfwi.h"

static int nm, nf, ww, rv;
static int *nn, **kk;
static float omega;
static sf_complex *tempx, *tempr, **ubac;

static float pi=3.141592653;

void fwi_init(int nt /* model size*/,
	      int *n /* model dimension */,
	      int ns /* number of shots */,
	      int shift /* shift */,
	      int **mask /* receiver mask */,
	      int water /* water layer depth */)
/*< initialize >*/
{
    /* model */
    nm = nt;
    nn = n;

    /* shot */
    nf = ns;

    /* shift */
    rv = shift;

    /* mask */
    kk = mask;

    /* water */
    ww = water;

    /* temp */
    tempx = sf_complexalloc(nt);
    tempr = sf_complexalloc(nt);

    /* background wavefield stored for re-usage */
    ubac = sf_complexalloc2(nt,ns);
}

void fwi_setup(float *sbac /* slowness background [n0*n1*n2] */,
	       float freq  /* frequency (Hz) */)
/*< setup preconditioner >*/
{
    /* helmholtz setup */
    /* NOTE: extra parameters: PML, alpha... */
    helm_setup(sbac,freq);

    /* angular frequency */
    omega = 2.*pi*freq;
}

void fwi_forward(sf_complex **source /* source */,
		 sf_complex **data /* data */,
		 sf_complex *dp /* misfit */)
/*< update data-misfit >*/
{
    int is, i0;

    /* loop over shot */
    for (is=0; is < nf; is++) {
	
	/* helmholtz solve */
	helm_solve(source[is],ubac[is]);

	/* receiver */
	for (i0=0; i0 < nn[1]*nn[2]; i0++) {
	    if (kk[is][i0])
		dp[i0+is*nn[1]*nn[2]] = data[is][rv+i0*nn[0]]-ubac[is][rv+i0*nn[0]];
	    else
		dp[i0+is*nn[1]*nn[2]] = sf_cmplx(0.,0.);
	}
    }
}

void fwi_operator(bool adj, bool add, int nx, int nr, sf_complex *x, sf_complex *r)
/*< fowrad/adjoint operator >*/
{
    int i0, iw, is;

    sf_cadjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dp solve ds */
	
	/* temp source */
	for (i0=0; i0 < nm; i0++) {
	    tempr[i0] = sf_cmplx(0.,0.);
	}
	
	/* loop over shot */
	for (is=0; is < nf; is++) {

	    /* receiver */
	    for (i0=0; i0 < nn[1]*nn[2]; i0++) {
		if (kk[is][i0])
		    tempr[rv+i0*nn[0]] = conjf(r[i0+is*nn[1]*nn[2]]);
		else
		    tempr[rv+i0*nn[0]] = sf_cmplx(0.,0.);
	    }

	    /* solve */
	    helm_solve(tempr,tempx);

	    for (i0=0; i0 < nm; i0++) {
		x[i0] -= omega*omega*conjf(tempx[i0]*ubac[is][i0]);
	    }
	}

	/* water layer */
	for (iw=0; iw < ww; iw++) {
	    for (i0=0; i0 < nn[1]*nn[2]; i0++) {
		x[iw+i0*nn[0]] = sf_cmplx(0.,0.);
	    }
	}
    } else {
	/* given ds solve dp */

	/* loop over shot */
	for (is=0; is < nf; is++) {

	    /* temp source */
	    for (i0=0; i0 < nm; i0++) {
		tempx[i0] = x[i0]*ubac[is][i0];
	    }

	    /* water layer */
	    for (iw=0; iw < ww; iw++) {
		for (i0=0; i0 < nn[1]*nn[2]; i0++) {
		    tempx[iw+i0*nn[0]] = sf_cmplx(0.,0.);
		}
	    }
	    
	    /* solve */
	    helm_solve(tempx,tempr);

	    /* receiver */
	    for (i0=0; i0 < nn[1]*nn[2]; i0++) {
		if (kk[is][i0])
		    r[i0+is*nn[1]*nn[2]] = -omega*omega*tempr[rv+i0*nn[0]];
		else
		    r[i0+is*nn[1]*nn[2]] = sf_cmplx(0.,0.);
	    }
	}
    }
}
