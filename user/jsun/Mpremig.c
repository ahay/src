/* Pseudo-spectral pre-stack source-receiver source independent diffraction imaging */
/*
   Copyright (C) 2009 University of Texas at Austin

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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psp.h"

int main(int argc, char* argv[])
{

    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs=1;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps, adj;
    float vref;

    bool roll; /* survey strategy */
    int offset,split;
    bool born;
    pspar par;
    int nx1, nz1; /*domain of interest*/
    float *vel,***dat,***dat_v,**wvfld1,**wvfld,*img,*imgs; /*velocity profile*/
    sf_file Fi,Fo,Fv,Fd_v,snaps; /* I/O files */
    sf_axis az,ax; /* cube axes */
    int shtbgn,shtend,shtnum,shtint;
    int ix,iz,is,which;
    bool justrec;

    bool diff;
    sf_file Fi1;
    float ***dat1;

    sf_init(argc,argv);

    if (!sf_getint("snap",&snap)) snap=0; /* interval for snapshots */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex fft */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if(!sf_getbool("abc",&abc)) abc=false; /* absorbing flag */
    if (!sf_getbool("roll", &roll)) roll=false; /*if n, receiver is independent of source location and gpl=nx*/
    if(!sf_getbool("born",&born)) born=false; /* born modeling flag */
    if(!sf_getbool("diff",&diff)) diff=false; /* diffraction imaging flag */
    if(!sf_getbool("justrec",&justrec)) justrec=false; /* just need full waveform record (no born or rtdm) */
    if (abc) {
        if(!sf_getint("nbt",&nbt)) sf_error("Need nbt!");
        if(!sf_getint("nbb",&nbb)) nbb = nbt;
        if(!sf_getint("nbl",&nbl)) nbl = nbt;
        if(!sf_getint("nbr",&nbr)) nbr = nbt;
        if(!sf_getfloat("ct",&ct)) sf_error("Need ct!");
        if(!sf_getfloat("cb",&cb)) cb = ct;
        if(!sf_getfloat("cl",&cl)) cl = ct;
        if(!sf_getfloat("cr",&cr)) cr = ct;
    } else {
        nbt = 0; nbb = 0; nbl = 0; nbr = 0;
        ct = 0; cb = 0; cl = 0; cr = 0;
    }
    if (!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getbool("ps",&ps)) ps=false; /* use pseudo-spectral */
    if (ps) sf_warning("Using pseudo-spectral...");
    else sf_warning("Using pseudo-analytical...");
    if (!sf_getbool("adj",&adj)) adj=false; /* use pseudo-spectral */
    if (justrec && adj) sf_error("No migration in justrec mode!");
    if (adj) sf_warning("RTM");
    else sf_warning("RTDM");
    if (!sf_getfloat("vref",&vref)) vref=1500; /* reference velocity (default using water) */

    /* setup I/O files */
    Fi = sf_input ("in");
    Fo = sf_output("out");
    Fv = sf_input("vel");
    if (adj) {
        gpl = -1;
        gpl_v = -1;
        sf_histint(Fi,"n1",&nt);
        sf_histfloat(Fi,"d1",&dt);
        sf_histint(Fi,"n2",&gpl);
        if (NULL!=sf_getstring("dat_v")) {
            Fd_v = sf_input("dat_v");
            sf_histint(Fd_v,"n2",&gpl_v);
        } else Fd_v = NULL;
        if (diff) Fi1 = sf_input("dat_2"); 
        else Fi1 = NULL;  
    } else {
        if (!sf_getint("nt",&nt)) sf_error("Need nt!");
        if (!sf_getfloat("dt",&dt)) sf_error("Need dt!");
        if (!sf_getint("gpl",&gpl)) gpl = -1; /* geophone length */
        if (!sf_getint("gpl_v",&gpl_v)) gpl_v = -1; /* geophone height */
    }
    if (!sf_getint("src",&src)) src=0; /* source type */
    //if (!sf_getint("n_srcs",&n_srcs)) n_srcs=1; /* source type */
    spx = sf_intalloc(n_srcs);
    spz = sf_intalloc(n_srcs);
    f0  = sf_floatalloc(n_srcs);
    t0  = sf_floatalloc(n_srcs);
    A   = sf_floatalloc(n_srcs);
    //if (!sf_getints("spx",spx,n_srcs)) sf_error("Need spx!"); /* shot position x */
    if (!sf_getints("spz",spz,n_srcs)) sf_error("Need spz!"); /* shot position z */
    if (!diff) {
        if (!sf_getfloats("f0",f0,n_srcs)) sf_error("Need f0! (e.g. 30Hz)");   /*  wavelet peak freq */
        if (!sf_getfloats("t0",t0,n_srcs)) sf_error("Need t0! (e.g. 0.04s)");  /*  wavelet time lag */
        if (!sf_getfloats("A",A,n_srcs)) sf_error("Need A! (e.g. 1)");     /*  wavelet amplitude */
    }
    if (!sf_getint("shtbgn", &shtbgn)) sf_error("Need shot starting location on grid!");
    if (!sf_getint("shtend", &shtend)) sf_error("Need shot ending location on grid!");
    if (!sf_getint("shtint", &shtint)) sf_error("Need shot interval on grid!");
    shtnum = (int)((shtend-shtbgn)/shtint) + 1;
    if (!sf_getint("which", &which)) which = 0;

    if (!sf_getint("gpx",&gpx)) gpx = -1; /* geophone position x */
    if (!sf_getint("gpz",&gpz)) gpz = -1; /* geophone position z */
    if (!sf_getint("gpx_v",&gpx_v)) gpx_v = -1; /* geophone position x */
    if (!sf_getint("gpz_v",&gpz_v)) gpz_v = -1; /* geophone position z */
    if (!sf_getint("offset",&offset)) offset = 0; /* nearest offset */
    if (!sf_getint("split",&split)) split = 1; /* receiver split */

    if (SF_FLOAT != sf_gettype(Fv)) sf_error("Need float input");

    /* Read/Write axes */
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);
    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;
    if (gpx==-1) gpx = nbl;
    if (gpz==-1) gpz = nbt;
    if (gpl==-1) gpl = nx1;
    if (gpx_v==-1) gpx_v = nbl;
    if (gpz_v==-1) gpz_v = nbt;
    if (gpl_v==-1) gpl_v = nz1;
    if (roll && gpl==nx1) sf_error("Rolling discrepency");

    if (adj) { /*output image*/
        sf_setn(az,nz1);
        sf_setn(ax,nx1);
        sf_oaxa(Fo,az,1);
        sf_oaxa(Fo,ax,2);
        sf_putint(Fo,"n3",1);
        sf_putfloat(Fo,"d3",shtint*dx);
        sf_putfloat(Fo,"o3",0.);
        sf_settype(Fo,SF_FLOAT);
    } else { /*output data*/
        sf_setn(ax,gpl);
        /*output horizontal data is mandatory*/
        sf_putint(Fo,"n1",nt);
        sf_putfloat(Fo,"d1",dt);
        sf_putfloat(Fo,"o1",0.);
        sf_putstring(Fo,"label1","Time");
        sf_putstring(Fo,"unit1","s");
        sf_oaxa(Fo,ax,2);
        sf_putint(Fo,"n3",shtnum);
        sf_putfloat(Fo,"d3",shtint*dx);
        sf_putfloat(Fo,"o3",0.);
        sf_putstring(Fo,"label3","Shot");
        sf_settype(Fo,SF_FLOAT);
        /*output vertical data is optional*/
        if (NULL!=sf_getstring("dat_v")) {
            Fd_v = sf_output("dat_v");
            sf_setn(ax,gpl_v);
            /*output horizontal data is mandatory*/
            sf_putint(Fd_v,"n1",nt);
            sf_putfloat(Fd_v,"d1",dt);
            sf_putfloat(Fd_v,"o1",0.);
            sf_putstring(Fd_v,"label1","Time");
            sf_putstring(Fd_v,"unit1","s");
            sf_oaxa(Fd_v,ax,2);
            sf_putint(Fd_v,"n3",shtnum);
            sf_putfloat(Fd_v,"d3",shtint*dx);
            sf_putfloat(Fd_v,"o3",0.);
            sf_putstring(Fd_v,"label3","Shot");
            sf_settype(Fd_v,SF_FLOAT);	
        } else Fd_v = NULL;
    }

    if (snap > 0) {
        snaps = sf_output("snaps");
        /* (optional) snapshot file */
        sf_setn(az,nz1);
        sf_setn(ax,nx1);
        sf_oaxa(snaps,az,1);
        sf_oaxa(snaps,ax,2);
        sf_putint(snaps,"n3",nt/snap);
        sf_putfloat(snaps,"d3",dt*snap);
        sf_putfloat(snaps,"o3",0.);
        sf_putstring(snaps,"label3","Time");
        sf_putstring(snaps,"unit3","s");
    } else snaps = NULL;

    par = (pspar) sf_alloc(1,sizeof(*par));
    vel = sf_floatalloc(nz*nx);
    dat = sf_floatalloc3(nt,gpl,shtnum);
    img = sf_floatalloc(nz1*nx1);
    if (adj) {
        imgs= sf_floatalloc(nz1*nx1);
        if (adj) {
            for (ix=0; ix<nx1; ix++)
                for (iz=0; iz<nz1; iz++)
                    imgs[ix*nz1+iz] = 0.;
        }
    } else imgs = NULL;
    if (NULL!=Fd_v) dat_v = sf_floatalloc3(nt,gpl_v,shtnum);
    else dat_v = NULL;
    if (snap>0) {
        wvfld1 = sf_floatalloc2(nx1*nz1,nt/snap);
        wvfld  = sf_floatalloc2(nx1*nz1,nt/snap);
    }
    else { wvfld1 = NULL; wvfld = NULL; }
    if (adj && diff) dat1 = sf_floatalloc3(nt,gpl,shtnum);
    else dat1 = NULL;

    sf_floatread(vel,nz*nx,Fv);
    if (adj) {
        sf_floatread(dat[0][0],gpl*nt*shtnum,Fi);
        if (NULL!=Fd_v) sf_floatread(dat_v[0][0],gpl_v*nt*shtnum,Fd_v);
        if (diff) sf_floatread(dat1[0][0],gpl*nt*shtnum,Fi1);
    } else {
        sf_floatread(img,nz1*nx1,Fi);
    }

    /*passing the parameters*/
    par->nx    = nx;  
    par->nz    = nz;
    par->dx    = dx;
    par->dz    = dz;
    par->n_srcs= n_srcs;
    par->spx   = spx;
    par->spz   = spz;
    par->gpz   = gpz;
    par->gpx   = gpx;
    par->gpl   = gpl;
    par->gpz_v = gpz_v;
    par->gpx_v = gpx_v;
    par->gpl_v = gpl_v;
    par->snap  = snap;
    par->cmplx = cmplx;
    par->pad1  = pad1;
    par->abc   = abc;
    par->nbt   = nbt;
    par->nbb   = nbb;
    par->nbl   = nbl;
    par->nbr   = nbr;
    par->ct    = ct;
    par->cb    = cb;
    par->cl    = cl;
    par->cr    = cr;
    par->src   = src;
    par->nt    = nt;
    par->dt    = dt;
    par->f0    = f0;
    par->t0    = t0;
    par->A     = A;
    par->verb  = verb;
    par->ps    = ps;
    par->vref  = vref;

    for (is=0; is<shtnum; is++){

        *spx = shtbgn + shtint*is;
        //par->spx = spx; //pointer
	if (roll) par->gpx = *spx+offset;

        sf_warning("Processing shot # %d/%d",is,shtnum-1);

        if (adj && diff) {

            sf_warning("Simultaneously propagating two receiver wavefields...");
            psp3(wvfld, wvfld1, dat[is], dat1[is], img, vel, par);


            for (ix=0; ix<nx1; ix++)
                for (iz=0; iz<nz1; iz++)
                    imgs[ix*nz1+iz] += img[ix*nz1+iz];

        } else {

            if (justrec) {

                if (NULL == dat_v)
                    psp(wvfld, dat[is], NULL, NULL, vel, par, false);
                else
                    psp(wvfld, dat[is], dat_v[is], NULL, vel, par, false);

            } else {

                sf_warning("Computing source wavefield ...");
                psp(wvfld1, NULL, NULL, NULL, vel, par, false);
                if (born) dt2v2(wvfld1, vel, par);
                sf_warning("Computing receiver wavefield ...");
                if (split>1) {
                    psp5(split, wvfld1, wvfld, dat[is], img, vel, par);
                } else {
                    if (NULL == dat_v)
                        psp2(wvfld1, wvfld, dat[is], NULL, img, vel, par, adj);
                    else
                        psp2(wvfld1, wvfld, dat[is], dat_v[is], img, vel, par, adj);
                }

                if (adj) {
                    for (ix=0; ix<nx1; ix++)
                        for (iz=0; iz<nz1; iz++)
                            imgs[ix*nz1+iz] += img[ix*nz1+iz];
                }

            }

        }

        if (snap>0 && is==which)
            sf_floatwrite(wvfld[0],nz1*nx1*nt/snap,snaps);
    }

    if (adj) {
        sf_floatwrite(imgs,nz1*nx1,Fo);
    } else {
        sf_floatwrite(dat[0][0],gpl*nt*shtnum,Fo);
        if (NULL!=Fd_v) sf_floatwrite(dat_v[0][0],gpl_v*nt*shtnum,Fd_v);
    }

    exit (0);
}
