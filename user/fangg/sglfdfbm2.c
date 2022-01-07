/* Forward and backward lowrank FD modeling on a staggered grid */
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
#include "pmlm.h"
#include "srcm.h"
#include "sglfdc.h"
#include "sglfdfbm2.h"
/*^*/
#ifndef _sglfdfbm2_h

typedef struct GeoPar {
    int   nx;
    int   nz;
    int   nxb;
    int   nzb;
    float dx;
    float dz;
    float ox;
    float oz;
    int   spx;
    int   spz;
    int   gp;
    int   gn;
    int   ginter;
    int   snpint;
   
} *geopar; /*geometry parameters*/
/*^*/

#endif

geopar creategeo(void)
/*< Create geometry used in RTM >*/
{
    geopar geop;
    geop = (geopar) sf_alloc(1, sizeof(*geop));
    return geop;
}

 
int sglfdfor2(float ***wavfld, float **rcd, bool verb,
              float **den, float **c11, 
              geopar geop, srcpar srcp, pmlpar pmlp)
/*< staggered grid lowrank FD forward modeling >*/
{
    /*caculate arrays*/
    float **txxn1, **txxn0, **vxn1, **vzn1, **vxn0, **vzn0;

    float **denx, **denz;
    /*grid index*/
    int nx, nz, nt, ix, iz, it;
    int nxb, nzb, snpint;
    int spx, spz, gp, gn, ginter;
    float dt;
    int pmlout, marg;
    bool freesurface;
    
    /* tmp variable */
    int wfit;

    nx = geop->nx;
    nz = geop->nz;
    nxb = geop->nxb;
    nzb = geop->nzb;

    spx = geop->spx;
    spz = geop->spz;
    gp  = geop->gp;
    gn  = geop->gn;
    ginter = geop->ginter;
    snpint = geop->snpint;

    nt = srcp->nt;
    dt = srcp->dt;

    pmlout = pmlp->pmlout;
    freesurface = pmlp->freesurface;
    marg   = getmarg();
    

    denx = sf_floatalloc2(nzb, nxb);
    denz = sf_floatalloc2(nzb, nxb);
    
    for (ix = 0; ix < nxb; ix++) {
	for ( iz= 0; iz < nzb; iz++) {
	    denx[ix][iz] = den[ix][iz];
	    denz[ix][iz] = den[ix][iz];
	}
    }
    /*den[ix+1/2][iz]*/
    for ( ix = 0; ix < nxb-1; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    denx[ix][iz] = (den[ix+1][iz] + den[ix][iz])*0.5;
	}
    }
    
    /*den[ix][iz+1/2]*/
    for ( ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb-1; iz++) {
	    denz[ix][iz] = (den[ix][iz+1] + den[ix][iz])*0.5;
	}
    }

    txxn1 = sf_floatalloc2(nzb, nxb);
    txxn0 = sf_floatalloc2(nzb, nxb);
    vxn1  = sf_floatalloc2(nzb, nxb);
    vzn1  = sf_floatalloc2(nzb, nxb);
    vxn0  = sf_floatalloc2(nzb, nxb);
    vzn0  = sf_floatalloc2(nzb, nxb);

    init_pml(nz, nx, dt, marg, pmlp);

    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    txxn1[ix][iz] = 0.0;
	 }
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    txxn0[ix][iz] = 0.0;
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vxn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vxn0[ix][iz] = 0.0;
	}
    } 
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vzn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vzn0[ix][iz] = 0.0;
	}
    }  
    for (it = 0; it < nt; it++) {
	for (ix = 0; ix < gn; ix++) {
	    rcd[ix][it] = 0.0;
	}
    }  

    /*Main loop*/
    wfit = 0;
    for (it = 0; it < nt; it++) {
	if (verb) sf_warning("it=%d/%d;", it, nt-1);
	
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++ ) {
	    for (iz = marg+pmlout; iz < nz+pmlout+marg; iz++) {
		vxn1[ix][iz] = vxn0[ix][iz] - dt/denx[ix][iz]*ldx(txxn0, ix, iz);
		vzn1[ix][iz] = vzn0[ix][iz] - dt/denz[ix][iz]*ldz(txxn0, ix, iz);
	    }
	}
	
	/*Velocity PML */
	pml_vxz(vxn1, vzn1, vxn0, vzn0, txxn0, denx, denz, ldx, ldz, freesurface);
	
	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
	    for ( iz = marg+pmlout; iz < nz+marg+pmlout; iz++) { 
		txxn1[ix][iz] = txxn0[ix][iz] - dt*c11[ix][iz]*(ldx(vxn1, ix-1, iz) + ldz(vzn1, ix, iz-1));
	    }
	}
	 
	/*Stress PML */
	pml_txx(txxn1, vxn1, vzn1, c11, ldx, ldz, freesurface);
	
	if ((it*dt)<=srcp->trunc) {
	    explsourcet(txxn1, srcp->wavelet, it, dt, spx+pmlout+marg, spz+pmlout+marg, nxb, nzb, srcp);
	}
	
	if ( it%snpint == 0 ) {
	    for ( ix = 0; ix < nx; ix++) 
		for ( iz = 0; iz<nz; iz++ )
		    wavfld[wfit][ix][iz] = txxn0[ix+pmlout+marg][iz+pmlout+marg];
	    wfit++;
	}
		 
	for ( ix =0 ; ix < gn; ix++) {
	    rcd[ix][it] = txxn0[ix*ginter+pmlout+marg][pmlout+marg+gp];
	    //sf_warning("rcd=%f ix=%d it=%d", rcd[ix][it], ix, it);
	}
	
	/*n1 -> n0*/
	time_step_exch(txxn0, txxn1, it);
	time_step_exch(vxn0, vxn1, it);
	time_step_exch(vzn0, vzn1, it);
	pml_tstep_exch(it);
	
    } /*Main loop*/
    sf_warning("wfit=%d",wfit);
    if (verb) sf_warning(".");
    return wfit;
    
}


int sglfdback2(float **img1, float **img2, float ***wavfld, float **rcd, 
               bool verb, float **den, float **c11, 
               geopar geop, srcpar srcp, pmlpar pmlp, sf_file Ftmpbwf)  
/*< staggered grid lowrank FD backward propagation + imaging >*/
{
    /*caculate arrays*/
    float **txxn1, **txxn0, **vxn1, **vzn1, **vxn0, **vzn0;
    
    float **denx, **denz;
    float **sill, **ccr;
    /*grid index*/
    int nx, nz, nt, ix, iz, it, gn, ginter;
    int nxb, nzb, snpint;
    int gp;
    float dt;
    int pmlout, marg;
    bool freesurface;

    /* tmp variable */
    int wfit;
    
    nx = geop->nx;
    nz = geop->nz;
    nxb = geop->nxb;
    nzb = geop->nzb;
    
    gp  = geop->gp;
    gn  = geop->gn;
    ginter = geop->ginter;
    snpint = geop->snpint;
    
    nt = srcp->nt;
    dt = srcp->dt;

    pmlout = pmlp->pmlout;
    freesurface = pmlp->freesurface;
    marg   = getmarg();

    denx = sf_floatalloc2(nzb, nxb);
    denz = sf_floatalloc2(nzb, nxb);

    sill = sf_floatalloc2(nz, nx);
    ccr  = sf_floatalloc2(nz, nx);
    
    for (ix = 0; ix < nxb; ix++) {
	for ( iz= 0; iz < nzb; iz++) {
	    denx[ix][iz] = den[ix][iz];
	    denz[ix][iz] = den[ix][iz];
	}
    }
    /*den[ix+1/2][iz]*/
    for ( ix = 0; ix < nxb-1; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    denx[ix][iz] = (den[ix+1][iz] + den[ix][iz])*0.5;
	}
    }
    /*den[ix][iz+1/2]*/
    for ( ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb-1; iz++) {
	    denz[ix][iz] = (den[ix][iz+1] + den[ix][iz])*0.5;
	}
    }
    
    txxn1 = sf_floatalloc2(nzb, nxb);
    txxn0 = sf_floatalloc2(nzb, nxb);
    vxn1  = sf_floatalloc2(nzb, nxb);
    vzn1  = sf_floatalloc2(nzb, nxb);
    vxn0  = sf_floatalloc2(nzb, nxb);
    vzn0  = sf_floatalloc2(nzb, nxb);

    init_pml(nz, nx, dt, marg, pmlp);
    
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    txxn1[ix][iz] = 0.0;
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    txxn0[ix][iz] = 0.0;
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vxn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vxn0[ix][iz] = 0.0;
	}
    } 
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vzn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nxb; ix++) {
	for (iz = 0; iz < nzb; iz++) {
	    vzn0[ix][iz] = 0.0;
	}
    }

    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    sill[ix][iz] = 0.0;
	}
    }
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    ccr[ix][iz] = 0.0;
	}
    }
        
    /*Main loop*/
    wfit = (int)(nt-1)/snpint;
    sf_warning("back wfit=%d",wfit);
    for (it = nt-1; it>=0; it--) {
	if  (verb) sf_warning("it=%d/%d;", it, nt-1);

	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
	    for ( iz = marg+pmlout; iz < nz+marg+pmlout; iz++) { 
		txxn0[ix][iz] = txxn1[ix][iz] + dt*c11[ix][iz]*(ldx(vxn1, ix-1, iz) + ldz(vzn1, ix, iz-1));
	    }
	}
	/*Stress PML */
	pml_txxb(txxn0, vxn1, vzn1, c11, ldx, ldz, freesurface);
	
	for (ix=0; ix<gn; ix++)  {
	    txxn0[ix*ginter+pmlout+marg][pmlout+marg+gp] = rcd[ix][it];
	}
	
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++ ) {
	    for (iz = marg+pmlout; iz < nz+pmlout+marg; iz++) {
		vxn0[ix][iz] = vxn1[ix][iz] + dt/denx[ix][iz]*ldx(txxn0, ix, iz);
		vzn0[ix][iz] = vzn1[ix][iz] + dt/denz[ix][iz]*ldz(txxn0, ix, iz);
	    }
	}

	/*Velocity PML */
	pml_vxzb(vxn1, vzn1, vxn0, vzn0, txxn0, denx, denz, ldx, ldz, freesurface);

	/*n1 -> n0*/
	time_step_exch(txxn1, txxn0, it);
	time_step_exch(vxn1, vxn0, it);
	time_step_exch(vzn1, vzn0, it);
	pml_tstep_exchb(it);
	
	if ( it%snpint == 0 ) {
	    for ( ix = 0; ix < nx; ix++) 
		sf_floatwrite(txxn0[ix+pmlout+marg]+pmlout+marg, nz, Ftmpbwf);
	}

	if (it%snpint == 0 ) {
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    ccr[ix][iz] += wavfld[wfit][ix][iz]*txxn0[ix+pmlout+marg][iz+pmlout+marg];
		}
	    }
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    sill[ix][iz] += wavfld[wfit][ix][iz]*wavfld[wfit][ix][iz];
		}
	    }
	    wfit--;
	}
    } /*Main loop*/

    
    /*crosscorrelation*/
    img1 = ccr; 
        
    /*crosscorrelation with source normalization*/
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    img2[ix][iz] = ccr[ix][iz]/(sill[ix][iz]+SF_EPS);
	}
    } 
    

	
    if (verb) sf_warning(".");
    return 0;

}


