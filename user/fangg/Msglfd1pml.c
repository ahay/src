/* 1-D Lowrank Finite-difference wave extrapolation */
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
    sf_file Fvel, Fden, Fsrc;
    sf_file Fwf/*wave field*/, Frec/*record*/, Fic/*Initial Condition*/; 
    sf_file FG, Fsx;
    /*I/O for MMS*/
    sf_file Fpsrc, Fvsrc, Fpint, Fvint;
    
    sf_axis at, ax;
    sf_axis icaxis;

    /*I/O arrays*/
    float *src; /*point source*/ 
    float *vel, *den, *c11, *record;
    float *ic;
    float *sxtmp;
    /*I/O for MMS*/
    float **psrc, **vsrc, *pint, *vint;
    
    /*Grid index variables*/
    int nx, nt, ix, it;
    int nxb;
    float dt, dx;

    /*caculate arrays*/
    float *txxn1, *txxn0, *vxn1, *vxn0;
    float *denx;

    /*source*/
    spara sp={0};
    bool  srcdecay, srcmms, inject;
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
    Fsrc = sf_input("in");  
    /*source wavelet*/
    Fvel    = sf_input("vel");
    /*velocity*/
    Fden    = sf_input("den");
    /*density*/
    Fwf     = sf_output("out");
    /*wavefield snap*/
    Frec    = sf_output("rec");
    /*record*/
    
    FG  = sf_input("G"); 
    Fsx = sf_input("sx");
    
    if (SF_FLOAT != sf_gettype(Fsrc)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fden)) sf_error("Need float input");
    
    /*parameters of source*/
    if (!sf_getbool("srcmms", &srcmms)) srcmms = false;
    /*source type: if y, use point source */
    if (!srcmms && !sf_getfloat("slx", &slx)) sf_error("Need slx input");
    /*source location in x */
    if (!srcmms && slx<0.0) sf_error("slx need >=0.0");
    if (!sf_getbool("srcdecay",&srcdecay)) srcdecay=false;
    /*source decay y=use*/
    if (!sf_getfloat("srctrunc",&srctrunc)) srctrunc=1000;
    /*source trunc time (s)*/
    if (!sf_getbool("inject", &inject)) inject = true;
    /* inject=y use inject source; inject=n use initial condition*/
    if (srcmms && !inject) sf_error("Initial condition and MMS are conflicted");
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
    if (!sf_histint(FG,"n2", &lenx)) sf_error("No n2= in input");

    /* Read/Write axes */   
    at = sf_iaxa(Fsrc,1); nt = sf_n(at); dt = sf_d(at);
    ax = sf_iaxa(Fvel,1); nxb = sf_n(ax); dx = sf_d(ax);   

    /*read FD coefficients*/
    G = sf_floatalloc2(nxb, lenx);
    sf_floatread(G[0], nxb*lenx, FG);

    /*read FD schemes*/
    sxtmp = sf_floatalloc(lenx);
    sx = sf_intalloc(lenx);
    sf_floatread(sxtmp, lenx, Fsx);
    mx = 0;
    for (ix=0; ix<lenx; ix++) {
	sx[ix] = (int)sxtmp[ix];
	mx = abs(sx[ix])>mx? abs(sx[ix]):mx;
    }
    marg = mx;
    
    nx = nxb - 2*pmlout - 2*marg;
        
    record = sf_floatalloc(nt);

    /*set wavefield axes*/     
    sf_setn(at, (int)(nt-1)/snapinter+1); /*set axis for snap file*/
    sf_setd(at,dt*snapinter);
    sf_setn(ax, nx);   
    sf_oaxa(Fwf,ax,1); 
    sf_oaxa(Fwf,at,2);

    /*set for record*/
    sf_setn(ax,1);
    
    /*read model*/
    vel  = sf_floatalloc(nxb);
    den  = sf_floatalloc(nxb);
    c11  = sf_floatalloc(nxb);
    denx = sf_floatalloc(nxb);
    
    sf_floatread(vel, nxb, Fvel);
    sf_floatread(den, nxb, Fden);
    for (ix = 0; ix < nxb; ix++) {
	c11[ix] = den[ix]*vel[ix]*vel[ix];
	denx[ix] = den[ix];
	if(c11[ix] <= 0.0) sf_warning("c11=%f: ix=%d ",c11[ix], ix);
    }
    /*den[ix+1/2]*/
    for ( ix = 0; ix < nxb-1; ix++) {
	denx[ix] = (den[ix+1] + den[ix])*0.5;
    }
    
    /*source and receiver location*/
    spx = (int)(slx/dx+0.5);
    gp  = (int)(gdep/dx+0.5);
    /*read source*/
    src = sf_floatalloc(nt);
    sf_floatread(src,nt,Fsrc);
       
    /*Initial Condition*/
    if (inject == false) {
	Fic = sf_input("ic");  
	/*initial condition*/
	if (SF_FLOAT != sf_gettype(Fic)) sf_error("Need float input of ic");
	icaxis = sf_iaxa(Fic, 1); 
	icnx = sf_n(icaxis);
	if (nx != icnx) sf_error("I.C. and velocity should be the same size.");
	ic = sf_floatalloc(nx);
	sf_floatread(ic, nx, Fic);	
    } else {
	ic = NULL;
    }

    /* Method of Manufactured Solution*/
    if (inject && srcmms) {
	Fpsrc = sf_input("presrc");
	Fvsrc = sf_input("velsrc");
	Fpint = sf_input("preinit");
	Fvint = sf_input("velinit");
		
	if (SF_FLOAT != sf_gettype(Fpsrc)) sf_error("Need float input of presrc");
	if (SF_FLOAT != sf_gettype(Fvsrc)) sf_error("Need float input of velsrc");
	if (SF_FLOAT != sf_gettype(Fpint)) sf_error("Need float input of preinit");
	if (SF_FLOAT != sf_gettype(Fvint)) sf_error("Need float input of velinit");
	
	psrc = sf_floatalloc2(nx, nt);
	vsrc = sf_floatalloc2(nx, nt);
	pint = sf_floatalloc(nx);
	vint = sf_floatalloc(nx);
	
	sf_floatread(psrc[0], nx*nt, Fpsrc);
	sf_floatread(vsrc[0], nx*nt, Fvsrc);
	sf_floatread(pint, nx, Fpint);
	sf_floatread(vint, nx, Fvint);
    } else {
	psrc = NULL;
	vsrc = NULL;
	pint = NULL;
	vint = NULL;
    }
            
    txxn1 = sf_floatalloc(nxb);
    txxn0 = sf_floatalloc(nxb);
    vxn1  = sf_floatalloc(nxb);
    vxn0  = sf_floatalloc(nxb);
 
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
    } 

    /* MAIN LOOP */
    it0 = 0;
    if (inject == false) {
	it0 = 1;
	for(ix = 0; ix < nx; ix++) {
	    txxn0[ix+marg+pmlout] = ic[ix];
	}
	sf_floatwrite(txxn0+pmlout+marg, nx, Fwf);
	record[0] = txxn0[pmlout+marg+gp];
    }
    
    /* MMS */
    if (inject && srcmms ) {
	it0 = 0;
	for (ix=0; ix <nx; ix++) {
	    txxn0[ix+marg+pmlout] = pint[ix]; /*P(x,0)*/
	    vxn0[ix+marg+pmlout]  = vint[ix]; /*U(x, -dt/2)*/
	}
    }
    
    for (it = it0; it < nt; it++) {
	if(verb) sf_warning("it=%d/%d;", it, nt-1);
	
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++ ) {
	    vxn1[ix] = vxn0[ix] - dt/denx[ix]*ldx(txxn0, ix);
	}

	/* MMS */
	if (inject && srcmms) 
	    for (ix = 0; ix < nx; ix++)
		vxn1[ix+marg+pmlout] += vsrc[it][ix]*dt;

	for (ix = 0; ix < nxb; ix++)
	    vxn0[ix] = vxn1[ix];

	/*Velocity PML */
	/* pml1_vxz(vxn1, vxn0, txxn0, denx, ldx, freesurface); */

	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
	    txxn1[ix] = txxn0[ix] - dt*c11[ix]*ldx(vxn1, ix-1);
	}
	
	/*Stress PML */
	/* pml1_txx(txxn1, vxn1, c11, ldx, freesurface); */

	if (inject) {
	    if (!srcmms && (it*dt)<=sp.trunc) {
		txxn1[marg+pmlout+spx] += src[it]*dt;
	    }
	    if (srcmms) {
		for (ix = 0; ix < nx; ix++) 
	    	    txxn1[ix+marg+pmlout] +=psrc[it][ix]*dt;
	    }
	}
	    
	if ( it%snapinter==0 ) {
	    sf_floatwrite(txxn0+pmlout+marg, nx, Fwf);
	}
	
	record[it] = txxn0[pmlout+marg+gp];
	for (ix=0; ix<nxb; ix++)
	    txxn0[ix] = txxn1[ix];

    }/*End of LOOP TIME*/
    
    sf_warning(".");
    sf_floatwrite(record, nt, Frec);
 
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglfd1pml is: %f seconds << ", duration);
    exit(0);
}    
    
    
    
    
    

    
