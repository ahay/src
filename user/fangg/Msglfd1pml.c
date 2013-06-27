/* 1-D Low Rank Finite-difference wave extrapolation */
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

static float **G;
static int *sx;
static int lenx;
static int marg;


static float ldx(float *data, int ix)
/*<Low rank finite difference : d/dx>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenx; il++) {
	res += 0.5*(-1*data[ix-sx[il]+1] + data[ix+sx[il]])*G[il][ix];
    }
    return res;
}

int main(int argc, char* argv[]) 
{
    clock_t tstart,tend;
    double  duration;
   
    /*flag*/
    bool verb;
    
    /*I/O*/
    sf_file fvel, fden, fsrc;
    sf_file fwf/*wave field*/, frec/*record*/, fic/*Initial Condition*/; 
    sf_file fG, fsx;
    
    sf_axis at, ax;
    sf_axis icaxis;

    /*I/O arrays*/
    float *srcp, **srcd; /*point source, distributed source*/ 
    float *vel, *den, *c11, *record;
    float *ic;
    float *sxtmp;
    
    /*Grid index variables*/
    int nx, nt, ix, it;
    int nxb;
    float dt, dx;

    /*caculate arrays*/
    float *txxn1, *txxn0, *vxn1, *vxn0;
    float *denx;

    /*source*/
    spara sp={0};
    bool  srcdecay, srcpoint, inject;
    float srctrunc;
    float slx;
    int   spx;

    /*PML*/
    int   pmlout, pmld0, decaybegin;
    int   decay;
    float gamma = GAMMA;
    int   mx; /*margin*/
    
    /*options*/
    bool  freesurface;
    int   gp;
    float gdep;
    int   snapinter;
      
    int it0;
    int icnx;
    
    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/

    /*Set I/O file*/
    fsrc = sf_input("in");  
    /*source wavelet*/
    fvel    = sf_input("vel");
    /*velocity*/
    fden    = sf_input("den");
    /*density*/
    fwf     = sf_output("out");
    /*wavefield snap*/
    frec    = sf_output("rec");
    /*record*/
    
    fG  = sf_input("G"); 
    fsx = sf_input("sx");
    
    if (SF_FLOAT != sf_gettype(fsrc)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");
    
    /*parameters of source*/
    if (!sf_getbool("srcpoint", &srcpoint)) srcpoint = true;
    /*source type: if y, use point source */
    if (srcpoint && !sf_getfloat("slx", &slx)) sf_error("Need slx input");
    /*source location in x */
    if (srcpoint && slx<0.0) sf_error("slx need >=0.0");
    if (!sf_getbool("srcdecay",&srcdecay)) srcdecay=false;
    /*source decay y=use*/
    if (!sf_getfloat("srctrunc",&srctrunc)) srctrunc=1000;
    /*source trunc time (s)*/
    if (!sf_getbool("inject", &inject)) inject = true;
    /* inject=y use inject source; inject=n use initial condition*/
    
     /*parameters of geometry*/
    if (!sf_getfloat("gdep", &gdep)) gdep = 0.0;
    /*depth of geophone */
    if (gdep <0.0) sf_error("gdep need to be >=0.0");
    /*source and receiver location*/

    if (!sf_getint("snapinter", &snapinter)) snapinter=1;
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
    if (!sf_histint(fG,"n2", &lenx)) sf_error("No n2= in input");

    /* Read/Write axes */   
    if (srcpoint) {
	at = sf_iaxa(fsrc,1);
    } else {
	at = sf_iaxa(fsrc,2);
    }
    nt = sf_n(at); dt = sf_d(at);
    ax = sf_iaxa(fvel,1); nxb = sf_n(ax); dx = sf_d(ax);   

    /*read FD coefficients*/
    G = sf_floatalloc2(nxb, lenx);
    sf_floatread(G[0], nxb*lenx, fG);

    /*read FD schemes*/
    sxtmp = sf_floatalloc(lenx);
    sx = sf_intalloc(lenx);
    sf_floatread(sxtmp, lenx, fsx);
    mx = 0;
    for (ix=0; ix<lenx; ix++) {
	sx[ix] = (int)sxtmp[ix];
	mx = abs(sx[ix])>mx? abs(sx[ix]):mx;
    }
    marg = mx;
    
    nx = nxb - 2*pmlout - 2*marg;
    
    /*read model*/
    vel  = sf_floatalloc(nxb);
    den  = sf_floatalloc(nxb);
    c11  = sf_floatalloc(nxb);
    denx = sf_floatalloc(nxb);
    
    sf_floatread(vel, nxb, fvel);
    sf_floatread(den, nxb, fden);
    for (ix = 0; ix < nxb; ix++) {
	c11[ix] = den[ix]*vel[ix]*vel[ix];
	denx[ix] = den[ix];
	if(c11[ix] <= 0.0) sf_warning("c11=%f: ix=%d ",c11[ix], ix);
    }
    /*den[ix+1/2]*/
    for ( ix = 0; ix < nxb-1; ix++) {
	denx[ix] = (den[ix+1] + den[ix])*0.5;
    }
    
    /*read source*/
    if (srcpoint) {
	srcp = sf_floatalloc(nt);
	sf_floatread(srcp,nt,fsrc);
    } else {
	srcd = sf_floatalloc2(nx,nt);
	sf_floatread(srcd[0], nx*nt, fsrc);
    }
        
    txxn1 = sf_floatalloc(nxb);
    txxn0 = sf_floatalloc(nxb);
    vxn1  = sf_floatalloc(nxb);
    vxn0  = sf_floatalloc(nxb);
        
    record = sf_floatalloc(nt);

    /*set wavefield axes*/     
    sf_setn(at, (int)(nt-1)/snapinter+1); /*set axis for snap file*/
    sf_setd(at,dt*snapinter);
    sf_setn(ax, nx);   
    sf_oaxa(fwf,ax,1); 
    sf_oaxa(fwf,at,2);

    /*set for record*/
    sf_oaxa(frec, at, 1);
    if (!srcpoint) {
	sf_setn(ax,1);
	sf_oaxa(frec, ax, 2);
    }

    spx = (int)(slx/dx+0.5);
    gp  = (int)(gdep/dx+0.5);
    
    /*Initial Condition*/
    if (inject == false) {
	fic = sf_input("ic");  
	/*initial condition*/
	if (SF_FLOAT != sf_gettype(fic)) sf_error("Need float input of ic");
	icaxis = sf_iaxa(fic, 1); 
	icnx = sf_n(icaxis);
	if (nx != icnx) sf_error("I.C. and velocity should be the same size.");
	ic = sf_floatalloc(nx);
	sf_floatread(ic, nx, fic);	
    }

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
    for (it = 0; it < nt; it++) {
	record[it] = 0.0;
    }  
    
    
    sp.trunc=srctrunc;
    sp.srange=10;
    sp.alpha=0.5;
    sp.decay=srcdecay?1:0;
   	
    if (verb) {
	sf_warning("============================");
	sf_warning("nx=%d nt=%d", nx, nt);
	sf_warning("dx=%f dt=%f", dx, dt);
	sf_warning("lenx=%d marg=%d pmlout=%d", lenx, marg, pmlout);
	sf_warning("srctrunc=%f srcdecay=%d", sp.trunc, sp.decay);
	sf_warning("slx=%f, spx=%d, gdep=%f gp=%d",slx,spx,gdep,gp);
	for(ix=0; ix<lenx; ix++){
	    sf_warning("[sxx]=[%d,] G=%f",sx[ix], G[ix][0]);
	}
	sf_warning("============================"); 
    } //End if

    /* MAIN LOOP */
    it0 = 0;
    if (inject == false) {
	it0 = 1;
	for(ix = 0; ix < nx; ix++) {
	    txxn0[ix+marg+pmlout] = ic[ix];
	    //vxn0[ix+marg+pmlout] = ic[ix];
	}
	sf_floatwrite(txxn0+pmlout+marg, nx, fwf);
	record[0] = txxn0[pmlout+marg+gp];
    }

    for (it = it0; it < nt; it++) {
	if(verb) sf_warning("it=%d/%d;", it, nt-1);
	
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++ ) {
	    vxn1[ix] = vxn0[ix] - dt/denx[ix]*ldx(txxn0, ix);
	}

	/*Velocity PML */
	pml1_vxz(vxn1, vxn0, txxn0, denx, ldx, freesurface);

	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
	    txxn1[ix] = txxn0[ix] - dt*c11[ix]*ldx(vxn1, ix-1);
	}
	
	/*Stress PML */
	pml1_txx(txxn1, vxn1, c11, ldx, freesurface);

	if (inject) {
	    if (srcpoint && (it*dt)<=sp.trunc) {
		//explsourcet1(txxn1, srcp, dt, it, spx+pmlout+marg, nt, &sp);
    		txxn1[marg+pmlout+spx] += srcp[it]*dt;
	    }
	    if (!srcpoint && (it*dt)<=sp.trunc) {
		for (ix = 0; ix < nx; ix++) 
	    	    txxn1[ix+marg+pmlout] +=srcd[it][ix]*dt;
	    }
	}
	    
	/*n1 -> n0*/
	time_step_exch1(txxn0, txxn1, it);
	time_step_exch1(vxn0, vxn1, it);
	pml1_tstep_exch(it);
	if ( it%snapinter==0 ) {
	    sf_floatwrite(txxn1+pmlout+marg, nx, fwf);
	}
	
	record[it] = txxn1[pmlout+marg+gp];
	
    }/*End of LOOP TIME*/
    
    sf_warning(".");
    sf_floatwrite(record, nt, frec);
 
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglfd1pml is: %f seconds << ", duration);
    exit(0);
}    
    
    
    
    
    

    
