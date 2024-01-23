/* PML bandary condtion for 1D staggered grid lowrank finite difference */
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
#include "pml1.h"

#define PMLOUT 30
#define PMLD0 300
#define DECAY_FLAG 1
#define DECAY_BEGIN 4000
#define GAMMA  0.1
/*^*/
/*variables from main function*/
static int pmlout, marg, nx;
static float dt;
static int pmld0, decay, decaybegin;
static float pmlgamma;
static float *pmldx, *decdx;
static float *txxn1x, *txxn0x;


void init_pml1(int vnx, float vdt,/* Modle size*/
	      int vpml,   int vmarg, /* pml boundary size and margine size */
	      int vpmld0, int vdecay, int vdecaybegin, float vgamma/* pml parameters*/)
/*<initialization>*/
{
    int nxb;
    int ix;
    
    pmlout = vpml;
    marg   = vmarg;
    nx  = vnx;
    dt  = vdt;
    
    pmld0 = vpmld0;
    decay = vdecay;
    decaybegin = vdecaybegin;
    pmlgamma = vgamma;
    
    nxb = nx+2*pmlout+2*marg;
    
    pmldx = sf_floatalloc(nxb);
    decdx = sf_floatalloc(nxb);
    
    txxn1x  =   sf_floatalloc(nxb);
    txxn0x  =   sf_floatalloc(nxb);
 
    /* PML and absorb */
    for (ix=0; ix<nxb; ix++){
	pmldx[ix] = 0.0;
	decdx[ix] = 1.0;
    }
    
    for(ix=marg;ix<marg+pmlout;ix++)
    {
	pmldx[ix]=pmld0*(marg+pmlout-ix)*(marg+pmlout-ix)/pmlout/pmlout;
	decdx[ix]=(-1)*pmlgamma*(marg+pmlout-ix)*(marg+pmlout-ix)/pmlout/pmlout+1;
    }
    
    for(ix=nx+marg+pmlout;ix<nx+marg+2*pmlout;ix++)
    {
	pmldx[ix]=pmld0*(ix-nx-marg-pmlout+1)*(ix-nx-marg-pmlout+1)/pmlout/pmlout;
	decdx[ix]=(-1)*pmlgamma*(ix-nx-marg-pmlout+1)*(ix-nx-marg-pmlout+1)/pmlout/pmlout+1;
    }
    
    for (ix = 0; ix <nxb; ix++){
	txxn1x[ix]  =   0.0;
	txxn0x[ix]  =   0.0;
    }
}


void pml1_vxz(float *vxn1, float *vxn0, 
	     float *txxn0, float *denx, 
	     float (*ldx)(float *, int), 
	     bool freesurface)
/*<velocity vx,vz  decay in pml>*/
{
    int ix;
    /*Velocity PML --top*/
    if (freesurface == false) {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    vxn1[ix]=((1-dt*pmldx[ix]/2)*vxn0[ix]-dt/denx[ix]*ldx(txxn0,ix))/(1+dt*pmldx[ix]/2);
	}
    } 
	
    /*Velocity PML  --bottom*/
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	vxn1[ix]=((1-dt*pmldx[ix]/2)*vxn0[ix]-dt/denx[ix]*ldx(txxn0,ix))/(1+dt*pmldx[ix]/2);
	
    }

 }

void pml1_txx(float *txxn1, float *vxn1, float *c11, 
	     float (*ldx)(float *, int),
             bool freesurface )
/*<stress decay in pml>*/
{
    int ix;
    /*Stress PML -- top*/
    if (freesurface == false) {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]-dt*c11[ix]*ldx(vxn1,ix-1))/(1+dt*pmldx[ix]/2);
	    txxn1[ix] = txxn1x[ix];
	}
    } else {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]-dt*0.0*ldx(vxn1,ix-1))/(1+dt*pmldx[ix]/2);
	    txxn1[ix]= txxn1x[ix];
	}
	
    }

    
    /*Stress PML -- bottom*/
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]-dt*c11[ix]*ldx(vxn1,ix-1))/(1+dt*pmldx[ix]/2);
	txxn1[ix] = txxn1x[ix];
    }	
}

void fdpml1_vxz(float *vxn1, float *vxn0, 
		float *txxn0, float *denx, 
		float dx, int oo, 
		float (*fdx)(float *, int, float, int),
		bool freesurface)
/*<velocity vx,vz  decay in pml --FD method>*/
{
    int ix;
    /*Velocity PML --top*/
    if (freesurface == false) {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    vxn1[ix] = ((1-dt*pmldx[ix]/2)*vxn0[ix] + dt/denx[ix]*fdx(txxn0, ix-1, dx, oo))/(1+dt*pmldx[ix]/2);
	    
    }
    } 
	
    /*Velocity PML  --bottom*/
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {	
	vxn1[ix] = ((1-dt*pmldx[ix]/2)*vxn0[ix] + dt/denx[ix]*fdx(txxn0, ix-1, dx, oo))/(1+dt*pmldx[ix]/2);
    }

}

void fdpml1_txx(float *txxn1, float *vxn1, float *c11, float dx, int oo,
		float (*fdx)(float *, int, float, int),
		bool freesurface )
/*<stress decay in pml>*/
{
    int ix;
    /*Stress PML -- top*/
    if (freesurface == false) {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]+dt*c11[ix]*fdx(vxn1,ix, dx, oo))/(1+dt*pmldx[ix]/2);
	    txxn1[ix] = txxn1x[ix];
	}
    } else {
	for (ix=marg; ix<marg+pmlout; ix++) {
	    txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]+dt*0.0*fdx(vxn1,ix,dx,oo))/(1+dt*pmldx[ix]/2);
	    txxn1[ix]= txxn1x[ix];
	}
    }

    
    /*Stress PML -- bottom*/
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	txxn1x[ix]=((1-dt*pmldx[ix]/2)*txxn0x[ix]+dt*c11[ix]*fdx(vxn1,ix,dx,oo))/(1+dt*pmldx[ix]/2);
	txxn1[ix] = txxn1x[ix];
    }	
}

void time_step_exch1(float *dn0, float *dn1, int it)
/*<exchange >*/
{
    int ix;
    for (ix = marg; ix<nx+pmlout*2+marg; ix++) {
	dn0[ix]  = dn1[ix];
    }
}
    

void pml1_tstep_exch(int it)
/*<pml exhange>*/
{
    time_step_exch1(txxn0x, txxn1x, it);
}

void pml1_close(void)
/*<free memory allocation>*/
{
    free(pmldx);
    free(decdx);
    free(txxn0x);
    free(txxn1x);
}




