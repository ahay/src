/* PML bandary condtion for staggered grid lowrank finite difference */
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
#include "pmlomp.h"

#define PMLOUT 30
#define PMLD0 300
#define DECAY_FLAG 1
#define DECAY_BEGIN 4000
#define GAMMA  0.1
/*^*/


/*variables from main function*/
static int pmlout, marg, nz, nx;
static float dt;
static int pmld0, decay, decaybegin;
static float pmlgamma;
static float *pmldx, *pmldz, *decdx, *decdz;
static float **txxn1x, **txxn1z, **txxn0x, **txxn0z;


void init_pml(int vnz,    int vnx, float vdt,/* Modle size*/
	      int vpml,   int vmarg, /* pml boundary size and margine size */
	      int vpmld0, int vdecay, int vdecaybegin, float vgamma/* pml parameters*/)
/*<initialization>*/
{
    int nxb, nzb;
    int ix, iz;
    
    pmlout = vpml;
    marg   = vmarg;
    nz  = vnz;
    nx  = vnx;
    dt  = vdt;
    
    pmld0 = vpmld0;
    decay = vdecay;
    decaybegin = vdecaybegin;
    pmlgamma = vgamma;
    
    nxb = nx+2*pmlout+2*marg;
    nzb = nz+2*pmlout+2*marg;

    pmldx = sf_floatalloc(nxb);
    pmldz = sf_floatalloc(nzb);
    decdx = sf_floatalloc(nxb);
    decdz = sf_floatalloc(nzb);

    txxn1x  =   sf_floatalloc2(nzb, nxb);
    txxn1z  =   sf_floatalloc2(nzb, nxb);
    txxn0x  =   sf_floatalloc2(nzb, nxb);
    txxn0z  =   sf_floatalloc2(nzb, nxb);
    

    /* PML and absorb */
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
    for (ix=0; ix<nxb; ix++){
	pmldx[ix] = 0.0;
	decdx[ix] = 1.0;
    }
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz<nzb; iz++){
	pmldz[iz] = 0.0;
	decdz[iz] = 1.0;
    }

#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
    for(ix=marg;ix<marg+pmlout;ix++)
    {
	pmldx[ix]=pmld0*(marg+pmlout-ix)*(marg+pmlout-ix)/pmlout/pmlout;
	decdx[ix]=(-1)*pmlgamma*(marg+pmlout-ix)*(marg+pmlout-ix)/pmlout/pmlout+1;
    }
    
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
    for(ix=nx+marg+pmlout;ix<nx+marg+2*pmlout;ix++)
    {
	pmldx[ix]=pmld0*(ix-nx-marg-pmlout+1)*(ix-nx-marg-pmlout+1)/pmlout/pmlout;
	decdx[ix]=(-1)*pmlgamma*(ix-nx-marg-pmlout+1)*(ix-nx-marg-pmlout+1)/pmlout/pmlout+1;
    }
    
    /* for pml_dz and damping dec_dz */
    for(iz=marg;iz<marg+pmlout;iz++)
    {
	pmldz[iz]=pmld0*(marg+pmlout-iz)*(marg+pmlout-iz)/pmlout/pmlout;
        decdz[iz]=(-1)*pmlgamma*(marg+pmlout-iz)*(marg+pmlout-iz)/pmlout/pmlout+1;
    }
    for(iz=nz+marg+pmlout;iz<nz+marg+2*pmlout;iz++)
    {
	pmldz[iz]=pmld0*(iz-nz-marg-pmlout+1)*(iz-nz-marg-pmlout+1)/pmlout/pmlout;
	decdz[iz]=(-1)*pmlgamma*(iz-nz-marg-pmlout+1)*(iz-nz-marg-pmlout+1)/pmlout/pmlout+1;
    }


#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix = 0; ix <nxb; ix++){
	for (iz = 0; iz<nzb; iz++) {
	    txxn1x[ix][iz]  =   0.0;
	    txxn1z[ix][iz]  =   0.0;
	    txxn0x[ix][iz]  =   0.0;
	    txxn0z[ix][iz]  =   0.0; 
	}
    }
}


void pml_vxz(float **vxn1, float **vzn1, float **vxn0, float **vzn0, 
	     float **txxn0, float **denx, float **denz,
	     float (*ldx)(float **, int, int), 
	     float (*ldz)(float **, int, int),
             bool freesurface)
/*<velocity vx,vz  decay in pml>*/
{
    int ix, iz;
    /*Velocity PML --top*/
    if (freesurface == false) {
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	    for (iz=marg; iz<marg+pmlout; iz++) {
		vxn1[ix][iz]=((1-dt*pmldx[ix]/2)*vxn0[ix][iz]-dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
		vzn1[ix][iz]=((1-dt*pmldz[iz]/2)*vzn0[ix][iz]-dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	    }
	}
    } 
	
    /*Velocity PML  --left*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<marg+pmlout; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    vxn1[ix][iz]=((1-dt*pmldx[ix]/2)*vxn0[ix][iz]-dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn1[ix][iz]=((1-dt*pmldz[iz]/2)*vzn0[ix][iz]-dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	}
    }
    /*Velocity PML  --right*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    vxn1[ix][iz]=((1-dt*pmldx[ix]/2)*vxn0[ix][iz]-dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn1[ix][iz]=((1-dt*pmldz[iz]/2)*vzn0[ix][iz]-dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	}
    }
    /*Velocity PML  --bottom*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout+nz; iz<nz+2*pmlout+marg; iz++) {
	    vxn1[ix][iz]=((1-dt*pmldx[ix]/2)*vxn0[ix][iz]-dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn1[ix][iz]=((1-dt*pmldz[iz]/2)*vzn0[ix][iz]-dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	}
    }
}

void pml_txx(float **txxn1, float **vxn1, float **vzn1, float **c11, 
	     float (*ldx)(float **, int, int), 
	     float (*ldz)(float **, int, int),
             bool freesurface )
/*<stress decay in pml>*/
{
    int ix, iz;
    /*Stress PML -- top*/
    if (freesurface == false) {
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	    for (iz=marg; iz<marg+pmlout; iz++) {
		txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
		txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
		txxn1[ix][iz] = txxn1x[ix][iz]+txxn1z[ix][iz];
	    }
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	    for (iz=marg; iz<marg+pmlout; iz++) {
		txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*0.0*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
		txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*0.0*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
		txxn1[ix][iz] = txxn1x[ix][iz]+txxn1z[ix][iz];
	    }
	}
    }

    /*Stress PML -- left*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<marg+pmlout; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn1[ix][iz] = txxn1x[ix][iz]+txxn1z[ix][iz];
	}
    }
    /*Stress PML -- right*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn1[ix][iz] = txxn1x[ix][iz]+txxn1z[ix][iz];
	}
    }
    /*Stress PML -- bottom*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout+nz; iz<nz+2*pmlout+marg; iz++) {
	    txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn1[ix][iz] = txxn1x[ix][iz]+txxn1z[ix][iz];
	}
    }	
}

void time_step_exch(float **dn0, float **dn1, int it)
/*<exchange >*/
{
    int ix, iz;
    if ( DECAY_FLAG == 1 && DECAY_BEGIN >= it){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix = marg; ix<nx+pmlout*2+marg; ix++) {
	    for (iz = marg; iz<nz+pmlout*2+marg; iz++) {
		dn0[ix][iz]  = dn1[ix][iz]*decdx[ix]*decdz[iz];
	    }
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix = marg; ix<nx+pmlout*2+marg; ix++) {
	    for (iz = marg; iz<nz+pmlout*2+marg; iz++) {
		dn0[ix][iz]  = dn1[ix][iz];
	    }
	}
    }
    
}

void pml_tstep_exch(int it)
/*<pml exhange>*/
{
    time_step_exch(txxn0x, txxn1x, it);
    time_step_exch(txxn0z, txxn1z, it);
}

void pml_close(void)
/*<free memory allocation>*/
{
    free(pmldx);
    free(pmldz);
    free(decdx);
    free(decdz);
    free(*txxn1z);
    free(txxn1z);
    free(*txxn0z);
    free(txxn0z);
    free(*txxn1x);
    free(txxn1x);
    free(*txxn0x);
    free(txxn1x);
}




