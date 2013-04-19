/* 1-D staggered Grid Finite-difference wave extrapolation */
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
#include "pml1.h"

static float fdx(float *xx, int ii, float ffdx, int oo)
{
    if(oo>=8)
	return ( 1.23409107e+000*(xx[ii+1] -xx[ii])-1.06649846e-001*(xx[ii+2]-xx[ii-1]) \
		 +2.30363667e-002*(xx[ii+3]-xx[ii-2])-5.34238560e-003*(xx[ii+4]-xx[ii-3]) \
		 +1.07727117e-003*(xx[ii+5] -xx[ii-4] )-1.66418878e-004*(xx[ii+6] -xx[ii-5] ) \
		 +1.70217111e-005*(xx[ii+7] -xx[ii-6] )-8.52346420e-007*(xx[ii+8] -xx[ii-7] ) )/ffdx;
    
    else if(oo==7)
	return ( 1.22860622e+000*(xx[ii+1] -xx[ii] )-1.02383852e-001*(xx[ii+2] -xx[ii-1] ) \
		 +2.04767704e-002*(xx[ii+3] -xx[ii-2] )-4.17893273e-003*(xx[ii+4] -xx[ii-3] ) \
		 +6.89453549e-004*(xx[ii+5] -xx[ii-4] )-7.69225034e-005*(xx[ii+6] -xx[ii-5] ) \
		 +4.23651475e-006*(xx[ii+7] -xx[ii-6] ) )/ffdx;
			
    else if(oo==6)
	return ( 1.22133636e+000*(xx[ii+1] -xx[ii] )-9.69314575e-002*(xx[ii+2] -xx[ii-1] ) \
		 +1.74476624e-002*(xx[ii+3] -xx[ii-2] )-2.96728952e-003*(xx[ii+4] -xx[ii-3] ) \
		 +3.59005398e-004*(xx[ii+5] -xx[ii-4] )-2.18478116e-005*(xx[ii+6] -xx[ii-5] ) )/ffdx;
    
    else if(oo==5)
	return ( 1.2112427*(xx[ii+1] -xx[ii] )-8.9721680e-2*(xx[ii+2] -xx[ii-1] )+1.3842773e-2*(xx[ii+3] -xx[ii-2] ) \
			-1.7656599e-3*(xx[ii+4] -xx[ii-3] )+1.1867947e-4*(xx[ii+5] -xx[ii-4] ) )/ffdx;
    
    else if(oo==4)
	return ( 1.19628906*(xx[ii+1] -xx[ii] )-7.97526042e-002*(xx[ii+2] -xx[ii-1] ) \
		 +9.57031250e-003*(xx[ii+3] -xx[ii-2] )-6.97544643e-004*(xx[ii+4] -xx[ii-3] ) )/ffdx;
    
    else if(oo==3)
	return ( 1.17187500*(xx[ii+1] -xx[ii] )-6.51041667e-002*(xx[ii+2] -xx[ii-1] ) \
		 +4.68750000e-003*(xx[ii+3] -xx[ii-2] ) )/ffdx;
    
    else if(oo==2)
	return ( 1.12500000*(xx[ii+1] -xx[ii] )-4.16666667e-002*(xx[ii+2] -xx[ii-1] ) )/ffdx;
    
    else if(oo==1)
	return (xx[ii+1] -xx[ii])/ffdx;
    else
    {sf_error("ERROE: \n");}
    
    return 0;
}



int main(int argc, char* argv[]) 
{
    clock_t tstart,tend;
    double duration;
    bool verb;
    int nx, nt, ix, it;
    int nxb;
    float dt, dx;
    float *txxn1, *txxn0, *vxn1, *vxn0;
    float *vel, *den, *c11, *source;
    float *denx, *denz;
    float *record;
    int spx, gdep;
    sf_file fvel, fden, fsource, fwf/*wave field*/, frec/*record*/;
    sf_axis at, ax;
    int marg;
    int snapinter;
    bool freesurface;
    spara sp={0};
    
    int pmlout, pmld0, decaybegin;
    int   decay;
    float gamma = GAMMA;
    
    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/

    /*Set I/O file*/
    fsource = sf_input("in");  /*source wavelet*/
    fvel    = sf_input("vel"); /*velocity*/
    fden    = sf_input("den"); /*density*/
    fwf     = sf_output("out");/*wavefield snap*/
     frec    = sf_output("rec"); /*record*/
    /* Read/Write axes */
    at = sf_iaxa(fsource, 1); nt = sf_n(at); dt = sf_d(at); 
    ax = sf_iaxa(fvel, 1); nxb = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(fwf, ax, 1);
    sf_oaxa(fwf, at, 2);

    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");
    
    if(!sf_getint("size",&marg)) marg=4;
    if (!sf_getint("gdep", &gdep)) gdep=0;
    /* recorder depth */
    if (!sf_getint("spx", &spx)) sf_error("Need spx input");
    /*source point in x */
    if (!sf_getint("snapinter", &snapinter)) snapinter=10;
    /* snap interval */
    if (!sf_getint("pmlsize", &pmlout)) pmlout=PMLOUT;
    /* size of PML layer */
    if (!sf_getint("pmld0", &pmld0)) pmld0=PMLD0;
    /* PML parameter */
    if (!sf_getint("decay",&decay)) decay=DECAY_FLAG;
    /* Flag of decay boundary condtion: 1 = use ; 0 = not use */
    if (!sf_getint("decaybegin",&decaybegin)) decaybegin=DECAY_BEGIN;
    /* Begin time of using decay boundary condition */
    if (!sf_getbool("freesurface", &freesurface)) freesurface=false;
    /*free surface*/
    
    nx = nxb - 2*pmlout - 2*marg;
    
    vel = sf_floatalloc(nxb);
    den = sf_floatalloc(nxb);
    c11 = sf_floatalloc(nxb);
    
    denx = sf_floatalloc(nx);
    denz = sf_floatalloc(nx);
    
    sf_floatread(vel, nxb, fvel);
    sf_floatread(den, nxb, fden);
    for (ix = 0; ix < nxb; ix++) {
	    c11[ix] = den[ix]*vel[ix]*vel[ix];
	    denx[ix] = den[ix];
	    if(c11[ix] == 0.0) sf_warning("c11=0: ix=%d iz%d", ix, iz);
	
    }
    /*den[ix+1/2][iz]*/
    for ( ix = 0; ix < nxb-1; ix++) {
	    denx[ix] = (den[ix+1] + den[ix])*0.5;
    }
    
    source = sf_floatalloc(nt);
    sf_floatread(source, nt, fsource);
    
    txxn1 = sf_floatalloc(nxb);
    txxn0 = sf_floatalloc(nxb);
    vxn1  = sf_floatalloc(nxb);
    vxn0  = sf_floatalloc(nxb);

    record = sf_floatalloc(nt);
    init_pml1(nx, dt, pmlout, marg, pmld0, decay, decaybegin, gamma);
    

    for (ix = 0; ix < nxb; ix++) {
	txxn1[ix] = 0.0;
    }
    for (ix = 0; ix < nxb; ix++) {
	txxn0[ix] = 0.0;
    }
    for (ix = 0; ix < nxb; ix++) {
	vxn1[ix] = 0.0;  
    }
    for (ix = 0; ix < nxb; ix++) {
	vxn0[ix] = 0.0;
    } 
    for (ix = 0; ix < nxb; ix++) {
	vzn1[ix] = 0.0;  
    }
    for (ix = 0; ix < nxb; ix++) {
	vzn0[ix] = 0.0;
    }  
    for (it = 0; it < nt; it++) {
	record[it] = 0.0;
    } 
   
    /* MAIN LOOP */
    sp.trunc=160;
    sp.srange=10;
    sp.alpha=0.5;
    sp.decay=1;
	
    sf_warning("============================");
    sf_warning("nx=%d  nt=%d", nx, nt);
    sf_warning("dx=%f  dt=%f", dx, dt);
    sf_warning("lenx=%d marg=%d pmlout=%d", lenx, marg, pmlout);
    
    
    for (it = 0; it < nt; it++) {
	sf_warning("it=%d;", it);
	if (it<=sp.trunc) {
	    explsourcet1(txxn0, source, it, spx+pmlout+marg, nxb, &sp);
	}
    
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++) {
		vxn1[ix] = vxn0[ix] + dt/den[ix][iz]*fdx(txxn0, ix-1, dx, oo);
	}
	
	/*Velocity PML*/
	pml1_vxz(vxn1, vxn0, txxn0, denx, ldx, freesurface);

	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
		txxn1[ix] = txxn0[ix] + dt*c11[ix][iz]*fdx(vxn1, ix, dx, oo);
	}
	
	/*Stress PML */
	pml1_txx(txxn1, vxn1, c11, ldx, freesurface);
	
	/*n1 -> n0*/
	time_step_exch1(txxn0, txxn1, it);
	time_step_exch1(vxn0, vxn1, it);
	pml1_tstep_exch(it);

	if ( it%snapinter==0 ) {
	    sf_floatwrite(txxn0+pmlout+marg, nx, fwf);
	}
	
	record[it] = txxn0[pmlout+marg+gdep];
	
    }/*End of LOOP TIME*/

    sf_warning(".");
    sf_floatwrite(record, nt, frec);
 
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsgfd1 is: %f seconds << ", duration);
    exit(0);
}    
    
    
    
    
    

    
