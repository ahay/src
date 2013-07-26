/* 2-D 4th-order Staggered Grid Finite-difference wave extrapolation */
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

static int M=4;

static float fdx(float **xx, int ii, int jj, float ffdx, int oo)
{
    if(oo>=4)
	return ( 1.19628906*(xx[ii+1][jj]-xx[ii][jj])-7.97526042e-002*(xx[ii+2][jj]-xx[ii-1][jj]) \
		 +9.57031250e-003*(xx[ii+3][jj]-xx[ii-2][jj])-6.97544643e-004*(xx[ii+4][jj]-xx[ii-3][jj]) )/ffdx;
    
    else if(oo==3)
	return ( 1.17187500*(xx[ii+1][jj]-xx[ii][jj])-6.51041667e-002*(xx[ii+2][jj]-xx[ii-1][jj]) \
		 +4.68750000e-003*(xx[ii+3][jj]-xx[ii-2][jj]) )/ffdx;
    
    else if(oo==2)
	return ( 1.12500000*(xx[ii+1][jj]-xx[ii][jj])-4.16666667e-002*(xx[ii+2][jj]-xx[ii-1][jj]) )/ffdx;
    
    else if(oo==1)
	return (xx[ii+1][jj]-xx[ii][jj])/ffdx;
    else
    {sf_error("ERROE: in fuction sfsglfd2_tfd!\n");}
    
    return 0;
}

static float fdz(float **xx, int ii, int jj, float ffdz, int oo)
{    
    if(oo>=4)
	return ( 1.19628906*(xx[ii][jj+1]-xx[ii][jj])-7.97526042e-002*(xx[ii][jj+2]-xx[ii][jj-1]) \
		 +9.57031250e-003*(xx[ii][jj+3]-xx[ii][jj-2])-6.97544643e-004*(xx[ii][jj+4]-xx[ii][jj-3]) )/ffdz;
    
    else if(oo==3)
	return ( 1.17187500*(xx[ii][jj+1]-xx[ii][jj])-6.51041667e-002*(xx[ii][jj+2]-xx[ii][jj-1]) \
		 +4.68750000e-003*(xx[ii][jj+3]-xx[ii][jj-2]) )/ffdz;
    
    else if(oo==2)
	return ( 1.12500000*(xx[ii][jj+1]-xx[ii][jj])-4.16666667e-002*(xx[ii][jj+2]-xx[ii][jj-1]) )/ffdz;
    
    else if(oo==1)
	return (xx[ii][jj+1]-xx[ii][jj])/ffdz;
    else
    {sf_error("ERROE: in fuction sfsglfd2_tfd!\n");}

    return 0;
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
    
    sf_file fvel, fden, fsource, fwf/*wave field*/; 
    sf_axis at, ax, az;
    int oo;

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

    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");
    
    if(!sf_getint("oo",&oo)) oo=4;
    if (!sf_getint("spx", &spx)) sf_error("Need spx input");
    /*source point in x */
    if (!sf_getint("spz", &spz)) sf_error("Need spz input");
    /* source point in z */
    
    vel = sf_floatalloc2(nz, nx);
    den = sf_floatalloc2(nz, nx);
    c11 = sf_floatalloc2(nz, nx);
    
    denx = sf_floatalloc2(nz, nx);
    denz = sf_floatalloc2(nz, nx);
    
    nxz = nx*nz; 
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
   
    /* MAIN LOOP */
    sp.trunc=0.2;
    sp.srange=10;
    sp.alpha=0.5;
    sp.decay=1.0;
	
    sf_warning("============================");
    sf_warning("nx=%d nz=%d nt=%d", nx, nz, nt);
    sf_warning("dx=%f dz=%f dt=%f", dx, dz, dt);

    
    
    for (it = 0; it < nt; it++) {
	sf_warning("it=%d;", it);
	if (it<=sp.trunc) {
	    explsourcet(txxn0, source, it, dt, spx, spz, nx, nz, &sp);
	}
    
	/*velocity*/
	for (ix = M; ix < nx-M; ix++ ) {
	    for (iz = M; iz < nz-M; iz++) {
		
		vxn1[ix][iz] = vxn0[ix][iz] + dt/den[ix][iz]*fdx(txxn0, ix-1, iz, dx, oo);
		vzn1[ix][iz] = vzn0[ix][iz] + dt/den[ix][iz]*fdz(txxn0, ix, iz, dz, oo);
	    }
	}
	
	/*Velocity PML*/

	/*Stress*/
	for (ix = M; ix < nx-M; ix++) {
	    for ( iz = M; iz < nz-M; iz++) { 
		txxn1[ix][iz] = txxn0[ix][iz] + dt*c11[ix][iz]*(fdx(vxn1, ix, iz, dx, oo) + fdz(vzn1, ix, iz-1, dz, oo));
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
    duration = (double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsgfd2 is: %f seconds << ", duration);
    exit(0);
}    
    
    
    
    
    

    
