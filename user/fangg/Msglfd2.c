/* 2-D Low Rank Finite-difference wave extrapolation */
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
#include <math.h>
#include <limits.h>
#include <time.h>
#include "source.h"


static float ***Gx, ***Gz;
static int *sxx, *sxz, *szx, *szz;
static int lenx, lenz;
static int M;


static float ldx(float **data, int ix, int iz)
/*<Low rank finite difference : d/dx>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenx; il++) {
	res += 0.5*(data[ix-sxx[il]][iz-sxz[il]] - data[ix+sxx[il]-1][iz+sxz[il]])*Gx[il][ix][iz];
    }
    return res;
}

static float ldz(float **data, int ix, int iz)
/*<Low rank finite difference : d/dz>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenz; il++) {
	res += 0.5*(data[ix-szx[il]][iz-szz[il]] - data[ix+szx[il]][iz+szz[il]-1])*Gz[il][ix][iz];
    }
    return res;
}

int main(int argc, char* argv[]) 
{
    clock_t tstart,tend;
    double duration;
    bool verb;
    int nx, nz, nxz, nt, ix, iz, it;
    float dt, dx, dz;
    float **txxn1, **txxn0, **vxn1, **vzn1, **vxn0, **vzn0;
    float **vel, **den, **c11, *source;
    float **denx, **denz;
    int spx, spz;
    
    int mx, mz; /*margin*/
 
    sf_file fvel, fden, fsource, fwf/*wave field*/; 
    sf_file fGx, fGz, fsxx, fsxz, fszx, fszz;
    float *sxxtmp, *sxztmp, *szxtmp, *szztmp;
    sf_axis at, ax, az;

    spara sp={0};
    
    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/

    /*Set I/O file*/
    fsource = sf_input("in");  /*source wavelet*/
    fvel    = sf_input("vel"); /*velocity*/
    fden    = sf_input("den"); /*density*/
    fwf     = sf_output("out");/*wavefield snap*/
    
    /* Read/Write axes */
    at = sf_iaxa(fsource, 1); nt = sf_n(at); dt = sf_d(at); 
    az = sf_iaxa(fvel, 1); nz = sf_n(az); dz = sf_d(az); 
    ax = sf_iaxa(fvel, 2); nx = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(fwf, az, 1);
    sf_oaxa(fwf, ax, 2);
    sf_oaxa(fwf, at, 3);

    fGx = sf_input("Gx");
    fGz = sf_input("Gz");
    fsxx = sf_input("sxx");
    fsxz = sf_input("sxz");
    fszx = sf_input("szx");
    fszz = sf_input("szz");
    
    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");
    
    if (!sf_getint("spx", &spx)) sf_error("Need spx input");
    /*source point in x */
    if (!sf_getint("spz", &spz)) sf_error("Need spz input");
    /* source point in z */
    if (!sf_histint(fGx, "n1", &nxz)) sf_error("No n1= in input");
    if (nxz != nx*nz) sf_error (" Need nxz = nx*nz");
    if (!sf_histint(fGx,"n2", &lenx)) sf_error("No n2= in input");
    if (!sf_histint(fGz,"n2", &lenz)) sf_error("No n2= in input");

    sxxtmp = sf_floatalloc(lenx);
    sxztmp = sf_floatalloc(lenx);
    szxtmp = sf_floatalloc(lenz);
    szztmp = sf_floatalloc(lenz);
    
    sxx = sf_intalloc(lenx);
    sxz = sf_intalloc(lenx);
    szx = sf_intalloc(lenz);
    szz = sf_intalloc(lenz);

    sf_floatread(sxxtmp, lenx, fsxx);
    sf_floatread(sxztmp, lenx, fsxz);
    sf_floatread(szxtmp, lenz, fszx);
    sf_floatread(szztmp, lenz, fszz);
    mx = 0; mz = 0;
    for (ix=0; ix<lenx; ix++) {
	sxx[ix] = (int)sxxtmp[ix];
	sxz[ix] = (int)sxztmp[ix];
	mx = abs(sxx[ix])>mx? abs(sxx[ix]):mx;
    }

    for (iz=0; iz<lenz; iz++) {
	szx[iz] = (int)szxtmp[iz];
	szz[iz] = (int)szztmp[iz];
	mz = abs(szz[iz])>mz? abs(szz[iz]):mz;
    }
    M = mx>mz?mx:mz;
    //M++;
    //M=9;

    Gx = sf_floatalloc3(nz, nx, lenx);
    Gz = sf_floatalloc3(nz, nx, lenz);
    sf_floatread(Gx[0][0], nz*nx*lenx, fGx);
    sf_floatread(Gz[0][0], nz*nx*lenz, fGz);
    free(sxxtmp); free(sxztmp); free(szxtmp); free(szztmp);
    vel = sf_floatalloc2(nz, nx);
    den = sf_floatalloc2(nz, nx);
    c11 = sf_floatalloc2(nz, nx);
    
    denx = sf_floatalloc2(nz, nx);
    denz = sf_floatalloc2(nz, nx);
    

    sf_floatread(vel[0], nxz, fvel);
    sf_floatread(den[0], nxz, fden);
    for (ix = 0; ix < nx; ix++) {
	for ( iz= 0; iz < nz; iz++) {
	    c11[ix][iz] = den[ix][iz]*vel[ix][iz]*vel[ix][iz];
	    denx[ix][iz] = den[ix][iz];
	    denz[ix][iz] = den[ix][iz];
	    if(c11[ix][iz] == 0.0) sf_warning("c11=0: ix=%d iz%d", ix, iz);
	}
    }
    /*den[ix+1/2][iz]*/
    for ( ix = 0; ix < nx-1; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    denx[ix][iz] = (den[ix+1][iz] + den[ix][iz])*0.5;
	}
    }
    
    /*den[ix][iz+1/2]*/
    for ( ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz-1; iz++) {
	    denz[ix][iz] = (den[ix][iz+1] + den[ix][iz])*0.5;
	}
    } 

    source = sf_floatalloc(nt);
    sf_floatread(source, nt, fsource);
    
    txxn1 = sf_floatalloc2(nz, nx);
    txxn0 = sf_floatalloc2(nz, nx);
    vxn1  = sf_floatalloc2(nz, nx);
    vzn1  = sf_floatalloc2(nz, nx);
    vxn0  = sf_floatalloc2(nz, nx);
    vzn0  = sf_floatalloc2(nz, nx);

    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    txxn1[ix][iz] = 0.0;
	 }
    }
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    txxn0[ix][iz] = 0.0;
	}
    }
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    vxn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    vxn0[ix][iz] = 0.0;
	}
    } 
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    vzn1[ix][iz] = 0.0;  
	}
    }
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    vzn0[ix][iz] = 0.0;
	}
    }  
   
    //sf_warning("I am here~~~~~~~~~~~~~~~~l205 ");
   
    /* MAIN LOOP */
    sp.trunc=160;
    sp.srange=10;
    sp.alpha=0.5;
    sp.decay=1;
	
    sf_warning("============================");
    sf_warning("nx=%d nz=%d nt=%d", nx, nz, nt);
    sf_warning("dx=%f dz=%f dt=%f", dx, dz, dt);
    sf_warning("lenx=%d lenz=%d M=%d", lenx, lenz, M);
    for(ix=0; ix<lenx; ix++){
	sf_warning("[sxx,sxz]=[%d,%d] Gx=%f",sxx[ix], sxz[ix], Gx[ix][0][0]);
    }
    for(ix=0; ix<lenz;ix++){
	sf_warning("[szx,szz]=[%d,%d] Gz=%f",szx[ix], szz[ix], Gz[ix][0][0]);
    } 
    
    for (it = 0; it < nt; it++) {
	if (it%20==0) sf_warning("it=%d;", it);
	if (it<=sp.trunc) {
	    explsourcet(txxn0, source, it, dt, spx, spz, nx, nz, &sp);
	}
    
	/*velocity*/
	for (ix = M; ix < nx-M; ix++ ) {
	    for (iz = M; iz < nz-M; iz++) {
		vxn1[ix][iz] = vxn0[ix][iz] - dt/denx[ix][iz]*ldx(txxn0, ix, iz);
		vzn1[ix][iz] = vzn0[ix][iz] - dt/denz[ix][iz]*ldz(txxn0, ix, iz);
	    }
	}

	/*Velocity PML*/

	/*Stress*/
	for (ix = M; ix < nx-M; ix++) {
	    for ( iz = M; iz < nz-M; iz++) { 
		txxn1[ix][iz] = txxn0[ix][iz] - dt*c11[ix][iz]*(ldx(vxn1, ix+1, iz) + ldz(vzn1, ix, iz+1));
	    }
	}
	
	for (ix = 0; ix<nx; ix++) {
	    for (iz = 0; iz < nz; iz++) {
		txxn0[ix][iz] = txxn1[ix][iz];
		vxn0[ix][iz] = vxn1[ix][iz];
		vzn0[ix][iz] = vzn1[ix][iz];
	    }
	}
	

	for ( ix = 0; ix < nx; ix++) {
	    sf_floatwrite(txxn0[ix], nz, fwf);
	}

	 
	
    }/*End of LOOP TIME*/
    sf_warning(".");

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglfd2 is: %f seconds << ", duration);

    exit(0);
}    
    
    
    
    
    

    
