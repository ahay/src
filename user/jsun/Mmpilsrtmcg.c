/* 2-D Low-rank One-step Least-square Pre-stack Reverse-Time-Migration using CG in the complex domain (both img and data are complex valued)
*/
/*
   Copyright (C) 2016 University of Texas at Austin

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
#include <mpi.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

/*automatically generated headers*/
#include "waveutils.h"

int psrtm(sf_complex*** record, sf_complex** imgsum, geopar geop);
void psrtm_lop(bool adj, bool add, int nm, int nd, sf_complex *mod, sf_complex *dat);
/*parameter structs*/
static geopar geop;

/*******************************************************/
/* main function */
int main(int argc, char* argv[]) 
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spz, gpz, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    sf_complex **ltf, **rtf;
    sf_complex **ltb, **rtb;
    sf_complex *ww;
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    /*misc*/
    int mode;
    int nzx, nx2, nz2, n2, nk;
    int ix, iz, it, is;
    int niter;

    /*Data/Image*/
    sf_complex ***record, **imginv;

    /*tmp*/
    int tmpint;

    /*I/O*/
    sf_file Fvel;
    sf_file left, right, leftb, rightb;
    sf_file Fsrc, Frcd/*source and record*/;
    sf_file Fimg;
    sf_file Fstart;

    /*axis*/
    sf_axis at, ax, az;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    sf_init(argc, argv);

    if(cpuid==0) sf_warning("numprocs=%d",numprocs);

    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/
    if (!sf_getint("niter",&niter)) niter=1;
    if (!sf_getint("mode", &mode)) mode = 0;
    if (!sf_getbool("illum", &illum)) illum=false; /*if n, no source illumination applied */
    if (!sf_getbool("roll", &roll)) roll=false; /*if n, receiver is independent of source location and gpl=nx*/
    /* source/receiver info */
    if (!sf_getint("shtbgn", &shtbgn)) sf_error("Need shot starting location on grid!");
    if (!sf_getint("sht0", &sht0)) sht0=shtbgn; /*actual shot origin on grid*/
    if (!sf_getint("shtend", &shtend)) sf_error("Need shot ending location on grid!");
    if (!sf_getint("shtint", &shtint)) sf_error("Need shot interval on grid!");
    shtnum = (int)((shtend-shtbgn)/shtint) + 1;
    shtnum0 = shtnum;
    if (shtnum%numprocs!=0) {
        shtnum += numprocs-shtnum%numprocs;
        if (verb) sf_warning("Total shot number is not divisible by total number of nodes! shunum padded from %d to %d.",shtnum0,shtnum); }
    if (!sf_getint("spz", &spz)) sf_error("Need source depth!");
    if (!sf_getint("gpz", &gpz)) sf_error("Need receiver depth!");
    if (roll) if (!sf_getint("gpl", &gpl)) sf_error("Need receiver length");
    if (!sf_getint("snapinter", &snpint)) snpint=1;     /* snap interval */
    /*--- parameters of source ---*/
    if (!sf_getfloat("srctrunc", &trunc)) trunc=0.4;
    if (!sf_getint("rectz", &rectz)) rectz=1;
    if (!sf_getint("rectx", &rectx)) rectx=1;
    if (!sf_getint("repeat", &repeat)) repeat=0;
    /* abc parameters */
    if (!sf_getint("top", &top)) top=40;
    if (!sf_getint("bot", &bot)) bot=40;
    if (!sf_getint("lft", &lft)) lft=40;
    if (!sf_getint("rht", &rht)) rht=40;

    /*Set I/O file*/
    Frcd = sf_input("--input"); /*record from elsewhere*/
    Fsrc  = sf_input("src");   /*source wavelet*/      
    Fimg  = sf_output("--output");
    left  = sf_input("left");
    right = sf_input("right");
    leftb  = sf_input("leftb");
    rightb = sf_input("rightb");
    Fvel  = sf_input("vel");  /*velocity - just for model dimension*/

    /*--- Axes parameters ---*/
    at = sf_iaxa(Fsrc, 1); nt = sf_n(at);  dt = sf_d(at);
    az = sf_iaxa(Fvel, 1); nzb = sf_n(az); dz = sf_d(az); oz = sf_o(az);
    ax = sf_iaxa(Fvel, 2); nxb = sf_n(ax); dx = sf_d(ax); ox = sf_o(ax);
    nzx = nzb*nxb;
    nz = nzb - top - bot;
    nx = nxb - lft - rht;
    if (!roll) gpl = nx; /* global survey setting */

    /* propagator matrices */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    nz2 = kiss_fft_next_fast_size(nzb*pad1);
    nx2 = kiss_fft_next_fast_size(nxb);
    nk = nz2*nx2; /*wavenumber*/
    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

    if (!sf_histint(leftb,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(leftb,"n2",&m2b))  sf_error("Need n2= in left");
    if (!sf_histint(rightb,"n1",&n2) || n2 != m2b) sf_error("Need n1=%d in right",m2b);
    if (!sf_histint(rightb,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

    /*check record data*/
    sf_histint(Frcd,"n1", &tmpint);
    if (tmpint != nt ) sf_error("Error parameter n1 in record!");
    sf_histint(Frcd,"n2", &tmpint);
    if (tmpint != gpl ) sf_error("Error parameter n2 in record!");
    sf_histint(Frcd,"n3", &tmpint);
    if (tmpint != shtnum0 ) sf_error("Error parameter n3 in record!");

    /*check starting model for cg*/
    if (NULL!=sf_getstring("start"))
        Fstart = sf_input("start");
    else
        Fstart = NULL;

    /*allocate memory*/
    ww=sf_complexalloc(nt);
    rr=sf_floatalloc(nzx);
    ltf = sf_complexalloc2(nzx,m2);
    rtf = sf_complexalloc2(m2,nk);
    ltb = sf_complexalloc2(nzx,m2b);
    rtb = sf_complexalloc2(m2b,nk);
    geop = (geopar) sf_alloc(1, sizeof(*geop));
    record = sf_complexalloc3(nt, gpl, shtnum);
    imginv = sf_complexalloc2(nz,nx);

    /*read from files*/
    sf_complexread(ww,nt,Fsrc);
    sf_complexread(ltf[0],nzx*m2,left);
    sf_complexread(rtf[0],m2*nk,right);
    sf_complexread(ltb[0],nzx*m2b,leftb);
    sf_complexread(rtb[0],m2b*nk,rightb);
    /*read data*/
    sf_complexread(record[0][0], shtnum0*gpl*nt, Frcd);
    if (shtnum0%numprocs!=0) {
#ifdef _OPENMP
#pragma omp parallel for private(is,ix,it)
#endif
        for (is=shtnum0; is<shtnum; is++)
            for (ix=0; ix<gpl; ix++)
                for (it=0; it<nt; it++)
                    record[is][ix][it] = sf_cmplx(0.,0.);
    }
    /*read starting model*/
    if (NULL != Fstart) {
        sf_complexread(imginv[0],nx*nz,Fstart);
    } else {
        /*transform the dimension to 1d*/
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
        for (iz=0; iz<nz; iz++) {
            for (ix=0; ix<nx; ix++) {
                imginv[ix][iz] = sf_cmplx(0.,0.);
            }
        }
    }

    /* output RSF files */
    if (cpuid==0) {
        sf_setn(az, nz);
        sf_setn(ax, nx);
        sf_oaxa(Fimg, az, 1);
        sf_oaxa(Fimg, ax, 2);
        sf_settype(Fimg,SF_COMPLEX);
    }

    /*close RSF files*/
    sf_fileclose(Fvel);
    sf_fileclose(Fsrc);
    sf_fileclose(left);
    sf_fileclose(right);
    sf_fileclose(leftb);
    sf_fileclose(rightb);

    /*load geopar elements*/
    geop->nx  = nx;
    geop->nz  = nz;
    geop->nxb = nxb;
    geop->nzb = nzb;
    geop->dx  = dx;
    geop->dz  = dz;
    geop->ox  = ox;
    geop->oz  = oz;
    geop->snpint = snpint;
    geop->spz = spz;
    geop->gpz = gpz;
    geop->gpl = gpl;
    geop->top = top;
    geop->bot = bot;
    geop->lft = lft;
    geop->rht = rht;
    geop->nt = nt;
    geop->dt = dt;
    geop->trunc = trunc;
    /*geop->adj   = adj; */
    geop->verb  = verb;
    geop->illum = illum;
    geop->m2 = m2;
    geop->m2b = m2b;
    geop->pad1 = pad1;
    /*pointers*/
    geop->ltf = ltf;
    geop->rtf = rtf;
    geop->ltb = ltb;
    geop->rtb = rtb;
    geop->ww  = ww;
    geop->rr  = rr;
    /*extra*/
    geop->roll = roll;
    geop->rectz=rectz;
    geop->rectx=rectx;
    geop->repeat=repeat;
    geop->sht0=sht0;
    geop->shtbgn=shtbgn;
    geop->shtend=shtend;
    geop->shtnum=shtnum;
    geop->shtnum0=shtnum0;
    geop->shtint=shtint;
    geop->cpuid=cpuid;
    geop->numprocs=numprocs;
    /*switch*/
    geop->mode=mode;

    sf_csolver(psrtm_lop,sf_ccgstep,nz*nx,nt*gpl*shtnum,imginv[0],record[0][0],niter,"verb",true,"end");
    //psrtm_lop(true,false,nz*nx,nt*gpl*shtnum,imginv[0],record[0][0]);

    if (cpuid==0) sf_complexwrite(imginv[0], nz*nx, Fimg);

    /*free memory*/
    free(geop);
    free(ww); free(rr);
    free(*ltf); free(ltf);
    free(*rtf); free(rtf);
    free(*ltb); free(ltb);
    free(*rtb); free(rtb);
    free(*imginv); free(imginv);
    free(**record); free(*record); free(record);

    MPI_Finalize();
    exit(0);
}

int psrtm(sf_complex*** record, sf_complex** imgsum, geopar geop)
    /*< low-rank one-step pre-stack RTM linear operator >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, lft; /*abc boundary*/
    int nt;
    float dt;
    bool adj; /* migration(adjoint) flag */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    /*pointers*/
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    sf_complex ***wvfld;
    sf_complex **tmprec, **img;
    float **sill;

    /*misc*/
    int ix, iz, is;
    int wfnt;
    float wfdt;

    clock_t tstart,tend;
    double duration;

    int shtcur;

    sf_complex *sendbuf, *recvbuf;

    /* passing variables */
    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    /*spx = geop->spx;*/
    spz = geop->spz;
    gpz = geop->gpz; 
    /*gpx = geop->gpx;*/
    gpl = geop->gpl;
    snpint = geop->snpint;
    top = geop->top; lft = geop->lft;
    nt = geop->nt;
    dt = geop->dt;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    rr = geop->rr;
    roll = geop->roll;
    rectz=geop->rectz;
    rectx=geop->rectx;
    repeat=geop->repeat;
    sht0=geop->sht0;
    shtbgn=geop->shtbgn;
    shtend=geop->shtend;
    shtnum=geop->shtnum;
    shtnum0=geop->shtnum0;
    shtint=geop->shtint;
    cpuid=geop->cpuid;
    numprocs=geop->numprocs;

    wfnt = (int)(nt-1)/snpint+1;
    wfdt = dt*snpint;

    /*allocate memory*/
    tmprec = sf_complexalloc2(nt, gpl);
    wvfld = sf_complexalloc3(nz, nx, wfnt);
    if (illum) sill = sf_floatalloc2(nz, nx);
    else sill = NULL;
    img = sf_complexalloc2(nz, nx);
    if (adj) { /*migration - initialize image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
        for (ix=0; ix<nx; ix++)
            for (iz=0; iz<nz; iz++)
                imgsum[ix][iz] = sf_cmplx(0.,0.);
    } else { /*modeling - initalize reflectivity*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
        for (ix=0; ix<nx; ix++)
            for (iz=0; iz<nz; iz++)
                img[ix][iz] = imgsum[ix][iz];
    }

    /* start the work */
    tstart = clock();

    for (is=0; is*numprocs<shtnum; is++){

        shtcur = is*numprocs+cpuid; // current shot index

        if (shtcur<shtnum0) {
            spx = shtbgn + shtint*(shtcur);
            if (roll)
                gpx = spx - (int)(gpl/2);
            else
                gpx = 0;
            geop->spx = spx;
            geop->gpx = gpx;

            if (verb) {
                sf_warning("============================");
                sf_warning("processing shot #%d", shtcur);
                sf_warning("nx=%d nz=%d nt=%d", geop->nx, geop->nz, geop->nt);
                sf_warning("nxb=%d nzb=%d ", geop->nxb, geop->nzb);
                sf_warning("dx=%f dz=%f dt=%f", geop->dx, geop->dz, geop->dt);
                sf_warning("top=%d bot=%d lft=%d rht=%d", geop->top, geop->bot, geop->lft, geop->rht);
                sf_warning("rectz=%d rectx=%d repeat=%d srctrunc=%f",rectz,rectx,repeat,geop->trunc);
                sf_warning("spz=%d spx=%d gpz=%d gpx=%d gpl=%d", spz, spx, gpz, gpx, gpl);
                sf_warning("snpint=%d wfdt=%f wfnt=%d ", snpint, wfdt, wfnt);
                sf_warning("sht0=%d shtbgn=%d shtend=%d shtnum0=%d shtnum=%d", sht0, shtbgn, shtend, shtnum0, shtnum);
                if (roll) sf_warning("Rolling survey!");
                else sf_warning("Global survey (gpl=nx)!");
                if (illum) sf_warning("Using source illumination!");
                else sf_warning("No source illumination!");
                sf_warning("============================");
            }

            /*generate reflectivity map*/
            reflgen(nzb, nxb, spz+top, spx+lft, rectz, rectx, repeat, rr);

            lrosfor2(wvfld, sill, tmprec, geop);
        }

        if (shtcur<shtnum0) {
            lrosback2(img, wvfld, sill, record[shtcur], geop);
            if (adj) { /*local image reduction*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
                for (ix=0; ix<nx; ix++) {
                    for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
                        imgsum[ix][iz] += img[ix][iz];
#else
                        imgsum[ix][iz] = sf_cadd(imgsum[ix][iz],img[ix][iz]);
#endif
                    }
                }
            }
        }

        if (!adj) {
            recvbuf = record[is*numprocs][0];
            sendbuf = tmprec[0];
            MPI_Allgather(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, MPI_COMM_WORLD); // record[is][ix][it] = tmprec[ix][it];
            /*Note that redundant shots will be received, but not outputed*/
        }

    } /*shot iteration*/

    MPI_Barrier(MPI_COMM_WORLD);

    /*image reduction*/
    if (adj) {
#if MPI_VERSION >= 2
        sendbuf = (sf_complex *) MPI_IN_PLACE;
#else /* will fail */
        sendbuf = NULL;
#endif 
        recvbuf = imgsum[0];
        MPI_Allreduce(sendbuf, recvbuf, nx*nz, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    if (verb)
        sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

    free(**wvfld); free(*wvfld); free(wvfld);
    free(*tmprec); free(tmprec);
    free(*img); free(img);
    if (illum) { free(*sill); free(sill);}

    return 0;
}

void psrtm_lop(bool adj, bool add, int nm, int nd, sf_complex *mod, sf_complex *dat)
    /*< lowrank onestep psrtm linear operator >*/
{
    sf_complex ***rec, **img;
    int nz, nx, nt, gpl, shtnum;
    int ix,is;

    nz = geop->nz;
    nx = geop->nx;
    nt = geop->nt;
    gpl = geop->gpl;
    shtnum = geop->shtnum;

    /*check the dimension*/
    if (nm!=nz*nx || nd!=nt*gpl*shtnum) sf_error("%s: wrong dimensions",__FILE__);
    sf_cadjnull(adj, add, nm, nd, mod, dat);

    /*allocate memory and transform the dimension*/
    img = (sf_complex**) sf_alloc (nx,sizeof(sf_complex*)); 
    img[0] = mod;
    for (ix=1; ix<nx; ix++) img[ix] = img[0]+ix*nz;

    rec = (sf_complex***) sf_alloc (shtnum,sizeof(sf_complex**));
    rec[0] = (sf_complex**) sf_alloc (gpl*shtnum,sizeof(sf_complex*));
    rec[0][0] = dat;
    for (ix=1; ix<gpl*shtnum; ix++) rec[0][ix] = rec[0][0]+ix*nt; 
    for (is=1; is<shtnum; is++) rec[is] = rec[0]+is*gpl;

    /* prestack migration */
    geop->adj=adj;
    psrtm(rec, img, geop);

    free(img);
    free(*rec); free(rec);
}
