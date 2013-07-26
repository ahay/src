/* PML bandary condtion for sglfd RTM */
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
#include "pmlm_OMP.h"

#ifndef _pmlm_OMP_h

#define PMLOUT 30
#define PMLD0 300
#define DECAY_FLAG true
#define DECAY_BEGIN 4000
#define GAMMA  0.1
/*^*/

typedef struct Pmlpar {
    int pmlout;
    int pmld0;
    bool decay;
    int decaybegin;
    float pmlgamma;
    bool freesurface;
} *pmlpar;
/*^*/

#endif

/*variables from main function*/
static int pmlout, marg, nz, nx;
static float dt;
static int pmld0, decaybegin;
static float pmlgamma;
static float *pmldx, *pmldz, *decdx, *decdz;
static float **txxn1x, **txxn1z, **txxn0x, **txxn0z;
static bool freesurface, decay;


pmlpar creatpmlpar(void)
/*< Create PML parameters >*/
{
    pmlpar pmlp;
    pmlp = (pmlpar) sf_alloc(1, sizeof(*pmlp));
    return pmlp;
}

void init_pml(int vnz,  int vnx, float vdt,/* Modle size*/
	      int vmarg, /* pml boundary size and margine size */
	      pmlpar pmlp)
/*< initialization >*/
{
    int nxb, nzb;
    int ix, iz;
    
    pmlout = pmlp->pmlout;
    marg   = vmarg;
    nz  = vnz;
    nx  = vnx;
    dt  = vdt;
    
    pmld0 = pmlp->pmld0;
    decay = pmlp->decay;
    decaybegin = pmlp->decaybegin;
    pmlgamma = pmlp->pmlgamma;
    freesurface = pmlp->freesurface;
    
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
    for (ix=0; ix<nxb; ix++){
	pmldx[ix] = 0.0;
	decdx[ix] = 1.0;
    }
    for (iz=0; iz<nzb; iz++){
	pmldz[iz] = 0.0;
	decdz[iz] = 1.0;
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
    
    //for pml_dz and damping dec_dz
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


void pml_vxzb(float **vxn1, float **vzn1, float **vxn0, float **vzn0, 
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
		vxn0[ix][iz]=((1-dt*pmldx[ix]/2)*vxn1[ix][iz]+dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
		vzn0[ix][iz]=((1-dt*pmldz[iz]/2)*vzn1[ix][iz]+dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	    }
	}
    } 
	
    /*Velocity PML  --left*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<marg+pmlout; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    vxn0[ix][iz]=((1-dt*pmldx[ix]/2)*vxn1[ix][iz]+dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn0[ix][iz]=((1-dt*pmldz[iz]/2)*vzn1[ix][iz]+dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	}
    }
    /*Velocity PML  --right*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    vxn0[ix][iz]=((1-dt*pmldx[ix]/2)*vxn1[ix][iz]+dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn0[ix][iz]=((1-dt*pmldz[iz]/2)*vzn1[ix][iz]+dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
	}
    }
    /*Velocity PML  --bottom*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout+nz; iz<nz+2*pmlout+marg; iz++) {
	    vxn0[ix][iz]=((1-dt*pmldx[ix]/2)*vxn1[ix][iz]+dt/denx[ix][iz]*ldx(txxn0,ix,iz))/(1+dt*pmldx[ix]/2);
	    vzn0[ix][iz]=((1-dt*pmldz[iz]/2)*vzn1[ix][iz]+dt/denz[ix][iz]*ldz(txxn0,ix,iz))/(1+dt*pmldz[iz]/2);
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

void pml_txxb(float **txxn0, float **vxn1, float **vzn1, float **c11, 
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
		txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
		txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
		txxn0[ix][iz] = txxn0x[ix][iz]+txxn0z[ix][iz];
	    }
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	    for (iz=marg; iz<marg+pmlout; iz++) {
		txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*0.0*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
		txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*0.0*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
		txxn0[ix][iz] = txxn0x[ix][iz]+txxn0z[ix][iz];
	    }
	}
    }

    /*Stress PML -- left*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<marg+pmlout; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn0[ix][iz] = txxn0x[ix][iz]+txxn0z[ix][iz];
	}
    }
    /*Stress PML -- right*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=nx+pmlout+marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout; iz<nz+pmlout+marg; iz++) {
	    txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn0[ix][iz] = txxn0x[ix][iz]+txxn0z[ix][iz];
	}
    }
    /*Stress PML -- bottom*/
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=marg; ix<nx+2*pmlout+marg; ix++) {
	for (iz=marg+pmlout+nz; iz<nz+2*pmlout+marg; iz++) {
	    txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*c11[ix][iz]*ldx(vxn1,ix-1,iz))/(1+dt*pmldx[ix]/2);
	    txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*c11[ix][iz]*ldz(vzn1,ix,iz-1))/(1+dt*pmldz[iz]/2);
	    txxn0[ix][iz] = txxn0x[ix][iz]+txxn0z[ix][iz];
	}
    }	
}



void time_step_exch(float **dn0, float **dn1, int it)
/*<exchange >*/
{
    int ix, iz;
    if ( decay && decaybegin >= it){
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

void pml_tstep_exchb(int it)
/*<pml exhange>*/
{
    time_step_exch(txxn1x, txxn0x, it);
    time_step_exch(txxn1z, txxn0z, it);
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




