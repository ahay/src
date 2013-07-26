/* 2-D Staggered Grid Lowrank Finite-difference RTM */
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
#include <time.h>

#include "sglfdfbm2_OMP.h"


void freertm()
{
    freesglfdcoe();
}

int main(int argc, char* argv[]) 
{
    clock_t tstart,tend;
    double duration;

    /*flag*/
    bool verb;
    
    /*I/O*/
    sf_file Ffvel, Ffden, Fbvel, Fbden;
    sf_file Fsrc,/*wave field*/ Frcd/*record*/;
    sf_file Ftmpwf, Ftmpbwf;
    sf_file Fimg;
    sf_file FGx, FGz, Fsxx, Fsxz, Fszx, Fszz;

    sf_axis at, ax, az;

    /*grid index variables*/
    int nx, nz, nxz, nt, wfnt;
    int ix, iz, it;
    int nxb, nzb;
    float dt, dx, dz, wfdt;
    float ox, oz;

    /*source/geophone location*/
    float slx, slz;
    int   spx, spz;
    float gdep;
    int   gp;

    /*SG LFD coefficient*/
    int lenx, lenz;
    int marg;

    /*Modle*/
    float **fvel, **fden, **fc11;
    float **bvel, **bden, **bc11;
    float ***wavefld, **record;
    float **img;
    int snpint;

    /*source*/
    bool srcdecay;
    int srcrange;
    float srctrunc;
    
    /*pml boundary*/
    int pmlout, pmld0, decaybegin;
    bool decay, freesurface;

    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/

    /*Set I/O file*/
    Fsrc   = sf_input("in");   /*source wavelet*/
    Ffvel  = sf_input("fvel"); /*forward velocity*/
    Ffden  = sf_input("fden"); /*forward density*/
    Fbvel  = sf_input("bvel"); /*backward velocity*/
    Fbden  = sf_input("bden"); /*backward velocity*/
    
    Ftmpwf = sf_output("tmpwf");/*wavefield snap*/
    Frcd   = sf_output("rec");   /*record*/
    Fimg   = sf_output("out");   /*Imaging*/
    Ftmpbwf = sf_output("tmpbwf");

    FGx = sf_input("Gx"); 
    FGz = sf_input("Gz");
    Fsxx = sf_input("sxx");
    Fsxz = sf_input("sxz");
    Fszx = sf_input("szx");
    Fszz = sf_input("szz");
    
    if (SF_FLOAT != sf_gettype(Ffvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Ffden)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fbvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fbden)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fsrc))  sf_error("Need float input");

   
    /*--- parameters of source ---*/
    srcpar srcp;
    srcp = createsrc();
    at = sf_iaxa(Fsrc, 1); nt = sf_n(at);  dt = sf_d(at);      
    if (!sf_getbool("srcdecay", &srcdecay)) srcdecay=SRCDECAY;
    /*source decay*/
    if (!sf_getint("srcrange", &srcrange)) srcrange=SRCRANGE;
    /*source decay range*/
    if (!sf_getfloat("srctrunc", &srctrunc)) srctrunc=SRCTRUNC;
    /*trunc source after srctrunc time (s)*/
    srcp->nt = nt; srcp->dt = dt; 
    srcp->decay = srcdecay; srcp->range=srcrange; srcp->trunc=srctrunc;
    loadsrc(srcp, Fsrc);

    /*--- parameters of SG LFD Coefficient ---*/
    ax = sf_iaxa(Ffvel, 2); nxb = sf_n(ax); dx = sf_d(ax); ox = sf_o(ax);
    az = sf_iaxa(Ffvel, 1); nzb = sf_n(az); dz = sf_d(az); oz = sf_o(az);

    if (!sf_histint(FGx, "n1", &nxz)) sf_error("No n1= in input");
    if (nxz != nxb*nzb) sf_error (" Need nxz = nxb*nzb");
    if (!sf_histint(FGx,"n2", &lenx)) sf_error("No n2= in input");
    if (!sf_histint(FGz,"n2", &lenz)) sf_error("No n2= in input");
    
    initsglfdcoe(nxb, nzb, lenx, lenz);
    loadcoe(nzb, nxb, FGx, FGz);
    loadschm(Fsxx, Fsxz, Fszx, Fszz);
    marg = getmarg();

    /* pml parameters */
    pmlpar pmlp;
    pmlp = creatpmlpar();
    if (!sf_getint("pmlsize", &pmlout)) pmlout=PMLOUT;
    /* size of PML layer */
    if (!sf_getint("pmld0", &pmld0)) pmld0=PMLD0;
    /* PML parameter */
    if (!sf_getbool("decay",&decay)) decay=DECAY_FLAG;
    /* Flag of decay boundary condtion: 1 = use ; 0 = not use */
    if (!sf_getint("decaybegin",&decaybegin)) decaybegin=DECAY_BEGIN;
    /* Begin time of using decay boundary condition */
    if (!sf_getbool("freesurface", &freesurface)) freesurface=false;
    /*free surface*/
    nx = nxb - 2*pmlout - 2*marg;
    nz = nzb - 2*pmlout - 2*marg;

    pmlp->pmlout =  pmlout;
    pmlp->pmld0 = pmld0;
    pmlp->decay = decay;
    pmlp->decaybegin = decaybegin;
    pmlp->freesurface = freesurface;

    /*Geometry parameters*/
    geopar geop;
    geop = creategeo();
 
    
    if (!sf_getfloat("slx", &slx)) sf_error("Need slx input");
    /*source location x */
    if (!sf_getfloat("slz", &slz)) sf_error("Need spz input");
    /* source location z */
    if (!sf_getfloat("gdep", &gdep)) gdep=0.0;
    /* recorder depth on grid*/
    if (gdep <0.0) sf_error("gdep need to be >=0.0");
    /*source and receiver location*/
    if (!sf_getint("snapinter", &snpint)) snpint=10;
    /* snap interval */
    spx = (int)((slx-ox)/dx+0.5);
    spz = (int)((slz-oz)/dz+0.5);
    gp  = (int)((gdep-oz)/dz+0.5);

    geop->nx  = nx;
    geop->nz  = nz;
    geop->nxb = nxb;
    geop->nzb = nzb;
    geop->dx  = dx;
    geop->dz  = dz;
    geop->ox  = ox;
    geop->oz  = oz;
    geop->snpint = snpint;
    geop->spx = spx;
    geop->spz = spz;
    geop->gp = gp;

    /* wavefield and record  */
    wfnt = (int)(nt-1)/snpint+1;
    wfdt = dt*snpint;
    record = sf_floatalloc2(nt, nx);
    wavefld = sf_floatalloc3(nz, nx, wfnt);
    
    /* read model */
    fvel = sf_floatalloc2(nzb, nxb);
    fden = sf_floatalloc2(nzb, nxb);
    fc11 = sf_floatalloc2(nzb, nxb);

    /*image*/
    img = sf_floatalloc2(nz, nx);

    
    sf_floatread(fvel[0], nxz, Ffvel);
    sf_floatread(fden[0], nxz, Ffden);
    for (ix = 0; ix < nxb; ix++) {
	for ( iz= 0; iz < nzb; iz++) {
	    fc11[ix][iz] = fden[ix][iz]*fvel[ix][iz]*fvel[ix][iz];
	    if(fc11[ix][iz] == 0.0) sf_warning("c11=0: ix=%d iz%d", ix, iz);
	}
    }

    bvel = sf_floatalloc2(nzb, nxb);
    bden = sf_floatalloc2(nzb, nxb);
    bc11 = sf_floatalloc2(nzb, nxb);
    
    sf_floatread(bvel[0], nxz, Fbvel);
    sf_floatread(bden[0], nxz, Fbden);
    for (ix = 0; ix < nxb; ix++) {
	for ( iz= 0; iz < nzb; iz++) {
	    bc11[ix][iz] = bden[ix][iz]*bvel[ix][iz]*bvel[ix][iz];
	    if(bc11[ix][iz] == 0.0) sf_warning("c11=0: ix=%d iz%d", ix, iz);
	}
    }
    if (verb) {
	sf_warning("============================");
	sf_warning("nx=%d nz=%d nt=%d", geop->nx, geop->nz, srcp->nt);
	sf_warning("dx=%f dz=%f dt=%f", geop->dx, geop->dz, srcp->dt);
	sf_warning("lenx=%d lenz=%d marg=%d pmlout=%d", lenx, lenz, marg, pmlout);
	sf_warning("srcdecay=%d srcrange=%d",srcp->decay,srcp->range);
	sf_warning("spx=%d spz=%d gp=%d snpint=%d", spx, spz, gp, snpint);
	sf_warning("wfdt=%f wfnt=%d ", wfdt, wfnt);
	sf_warning("============================");
    }

    /* write record */
    sf_setn(ax, nx);
    sf_setn(az, nz);
    
    sf_oaxa(Frcd, at, 1);
    sf_oaxa(Frcd, ax, 2);

    /*write temp wavefield */
    sf_setn(at, wfnt);
    sf_setd(at, wfdt);

    sf_oaxa(Ftmpwf, az, 1);
    sf_oaxa(Ftmpwf, ax, 2);
    sf_oaxa(Ftmpwf, at, 3);

    /*write temp wavefield */
    sf_oaxa(Ftmpbwf, az, 1);
    sf_oaxa(Ftmpbwf, ax, 2);
    sf_oaxa(Ftmpbwf, at, 3);

    /*write image*/
    sf_oaxa(Fimg, az, 1);
    sf_oaxa(Fimg, ax, 2);

    sglfdfor2(wavefld, record, verb, fden, fc11, geop, srcp, pmlp);
    sf_warning("Wavefield nt=%d", wfnt);
    sglfdback2(img, wavefld, record, verb, bden, bc11, geop, srcp, pmlp, Ftmpbwf);
   
    
    for (ix=0; ix<nx; ix++) 
	sf_floatwrite(record[ix], nt, Frcd);
    
    for (it=0; it<wfnt; it++)
	for ( ix=0; ix<nx; ix++)
	    sf_floatwrite(wavefld[it][ix], nz, Ftmpwf);

    for (ix=0; ix<nx; ix++) 
	sf_floatwrite(img[ix], nz, Fimg);
    
    
    freertm();
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglfd2 is: %f seconds << ", duration);
    exit(0);

}



