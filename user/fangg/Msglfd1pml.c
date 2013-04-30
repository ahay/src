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
	res += 0.5*(data[ix-sx[il]] - data[ix+sx[il]-1])*G[il][ix];
    }
    return res;
}

int main(int argc, char* argv[]) 
{
    clock_t tstart,tend;
    double duration;
    bool verb;
    int nx, nt, ix, it, it0;
    int nxb;
    float dt, dx;
    float *txxn1, *txxn0, *vxn1, *vxn0;
    float *vel, *den, *c11, *source;
    float *denx;
    float *record;
    bool freesurface;
    bool inject = false;
    int spx, gdep;
    int snapinter;
    int mx; /*margin*/
 
    sf_file fvel, fden, fsource, fwf/*wave field*/, frec/*record*/, fic/*Initial Condition*/; 
    sf_file fG, fsx;
    float *sxtmp;
    sf_axis at, ax;

    spara sp={0};
    bool srcdecay;
    int srctrunc;

    int pmlout, pmld0, decaybegin;
    int   decay;
    float gamma = GAMMA;

    int icnx;
    sf_axis icaxis;
    float *ic;

    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/

    /*Set I/O file*/
    fsource = sf_input("in");  
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
        
    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");
        
    if (!sf_getint("spx", &spx)) sf_error("Need spx input");
    /*source point in x */
    if (!sf_getint("gdep", &gdep)) gdep=0;
    /* recorder depth */
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
    if (!sf_histint(fG,"n2", &lenx)) sf_error("No n2= in input");
    if (!sf_getbool("srcdecay",&srcdecay)) srcdecay=true;
    /*source decay y=use*/
    if (!sf_getint("srctrunc",&srctrunc)) srctrunc=300;
    /*source trunc*/
    if (!sf_getbool("inject", &inject)) inject = true;
    /* inject=y use inject source; inject=n use initial condition*/
    
    /* Read/Write axes */
    at = sf_iaxa(fsource, 1); nt = sf_n(at); dt = sf_d(at); 
    ax = sf_iaxa(fvel,1); nxb = sf_n(ax); dx = sf_d(ax);
    
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
    
    /*set axis for record file*/
    //sf_setn(at, nt);
    sf_setn(ax, nx);
    sf_oaxa(frec, at, 1);
    

    /*set axis for snap file*/
    sf_setn(at, (int)(nt-1)/snapinter+1);
    sf_setd(at,dt*snapinter);
    sf_oaxa(fwf, ax, 1);
    sf_oaxa(fwf, at, 2);
    
    G = sf_floatalloc2(nxb, lenx);
    sf_floatread(G[0], nxb*lenx, fG);
    free(sxtmp);
    vel = sf_floatalloc(nxb);
    den = sf_floatalloc(nxb);
    c11 = sf_floatalloc(nxb);
    
    denx = sf_floatalloc(nxb);
    
    sf_floatread(vel, nxb, fvel);
    sf_floatread(den, nxb, fden);
    for (ix = 0; ix < nxb; ix++) {
	c11[ix] = den[ix]*vel[ix]*vel[ix];
	denx[ix] = den[ix];
	if(c11[ix] <= 0.0) sf_warning("c11=%f: ix=%d ",c11[ix], ix);
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
    
    
    /*Initial Condition*/
    if (inject == false) {
	fic     = sf_input("ic");  
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
   	
    sf_warning("============================");
    sf_warning("nx=%d nt=%d", nx, nt);
    sf_warning("dx=%f dt=%f", dx, dt);
    sf_warning("lenx=%d marg=%d pmlout=%d", lenx, marg, pmlout);
    sf_warning("srctrunc=%d srcdecay=%d", sp.trunc, sp.decay);
    for(ix=0; ix<lenx; ix++){
	sf_warning("[sxx]=[%d,] G=%f",sx[ix], G[ix][0]);
    }
    
    
    /* MAIN LOOP */
    for (it = 0; it < nt; it++) {
	sf_warning("it=%d/%d;", it, nt);
	if (inject ==true && it<=sp.trunc) {
	    explsourcet1(txxn0, source, it, spx+pmlout+marg, nxb, &sp);
	}
	
	/*velocity*/
	for (ix = marg+pmlout; ix < nx+pmlout+marg; ix++ ) {
	    vxn1[ix] = vxn0[ix] - dt/denx[ix]*ldx(txxn0, ix);
	}

	/*Velocity PML */
	pml1_vxz(vxn1, vxn0, txxn0, denx, ldx, freesurface);

	/*Stress*/
	for (ix = marg+pmlout; ix < nx+marg+pmlout; ix++) {
	    txxn1[ix] = txxn0[ix] - dt*c11[ix+1]*ldx(vxn1, ix+1);
	}
	
	/*Stress PML */
	pml1_txx(txxn1, vxn1, c11, ldx, freesurface);

	/*n1 -> n0*/
	time_step_exch1(txxn0, txxn1, it);
	time_step_exch1(vxn0, vxn1, it);
	pml1_tstep_exch(it);
	
	
	if (it ==0 && inject == false) {
	    for(ix = 0; ix < nx; ix++) {
		txxn0[ix+marg+pmlout] = ic[ix];
		//vxn0[ix+marg+pmlout] = ic[ix];
	    }
	}
	
       	if ( it%snapinter==0 ) {
	    sf_floatwrite(txxn0+pmlout+marg, nx, fwf);
	}
	
	record[it] = txxn0[pmlout+marg+gdep];
	
    }/*End of LOOP TIME*/
    sf_warning(".");
    sf_floatwrite(record, nt, frec);
 
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglfd1pml is: %f seconds << ", duration);
    exit(0);
}    
    
    
    
    
    

    
