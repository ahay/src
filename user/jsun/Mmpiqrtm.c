/* 2-D Low-rank One-step Least Pre-stack Reverse-Time-Migration in the complex domain (both img and data are complex valued)
     img :  crosscorrelation with source normalization (stdout)
*/
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
#include "compensate.h"

int psrtm(sf_complex*** record, sf_complex** imgsum, geopar geop);
int psrtm_mov(sf_complex*** record, sf_complex** imgsum, geopar geop, sf_file Ftmpwf, sf_file Ftmpwfb, int shtid);
int psqrtm_sbs(sf_complex*** record, sf_complex** imgsum, geopar geop, bool sdiv);
int psqrtm_com(sf_complex*** record, sf_complex** imgsum, geopar geop, bool freq_scal, sf_file Ftmpwf, sf_file Ftmpwfb, int shtid);
int psqrtm_dec(sf_complex*** record, sf_complex** imgsum, geopar geop);

/* for stable division */
static float eps, perc;
static int dims[2], rect[2], niter;
static float reg;

/*******************************************************/
/* main function */
int main(int argc, char* argv[]) 
{
    /*parameter structs*/
    geopar geop;

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
    bool adj; /* migration(adjoint) flag */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    sf_complex **ltf, **rtf;
    sf_complex **ltb, **rtb;
    sf_complex *ww;
    float *rr;
    /*extras*/
    bool roll,freq_scal; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    /*misc*/
    int mode;
    int nzx, nx2, nz2, n2, nk;
    int ix, iz, it, is;
    int wfnt;
    float wfdt;
    int stable;
    bool sdiv; /*smooth division by least squares*/
    int shtid; /*output wavefield corresponding shot id*/
    float vmax, veps; /* for stable division */

    /*Data/Image*/
    sf_complex ***record, **imgsum;
    sf_complex *img_visc, *img_disp;
    float *ratio;
    float *img_visc_f, *img_disp_f;

    /*tmp*/
    int tmpint;

    /*I/O*/
    sf_file Fvel;
    sf_file left, right, leftb, rightb;
    sf_file Fsrc, Frcd/*source and record*/;
    sf_file Fimg;
    sf_file Fnorm;
    sf_file Ftmpwf, Ftmpwfb;

    /*axis*/
    sf_axis at, ax, az, as;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    sf_init(argc, argv);

    if(cpuid==0) sf_warning("numprocs=%d",numprocs);

    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/
    if (!sf_getbool("adj", &adj)) adj=true; /*migration*/
    if (!sf_getint("mode", &mode)) {
      if (adj) mode = 1; /*select propagator*/
      else mode = 0;
    }
    if (verb && cpuid==0) sf_warning("mode=%d",mode);
    if (!sf_getbool("illum", &illum)) illum=false; /*if n, no source illumination applied */
    if (!sf_getbool("roll", &roll)) roll=false; /*if n, receiver is independent of source location and gpl=nx*/
    if (!sf_getint("stable", &stable)) stable=0; /*stable = 0 -> conventional imaging condition; 1 -> stable imaging condition for Q-compensation with global nomalization; 2 -> shot-by-shot normalization; 3 -> snapshot-by-snapshot compensation (most intensive); 4 -> deconvolution imaging condition */
    if (adj && stable!=0) {
      if (cpuid==0) {
	sf_warning("For stable imaging condition, mode parameter is preceded!");
	if (stable==1) sf_warning("global normalization!");
	else if (stable==2) sf_warning("shot by shot normalization!");
	else if (stable==3) sf_warning("snapshot-by-snapshot compensation, illum option is not recommended!");
	else if (stable==4) sf_warning("deconvolution imaging condition, illum parameter is disabled!");
	else sf_error("Please specify a correct stable parameter (0,1,2,3)!");
      }
      if (stable==1 || stable==2) {
          if (!sf_getbool("sdiv", &sdiv)) sdiv=false; /*smooth division*/
          if (sdiv) {
              if (!sf_getfloat("reg", &reg)) reg=0.0f; /*regularization*/
              if (!sf_getint("rect1",rect)) rect[0]=2;
              if (!sf_getint("rect2",rect+1)) rect[1]=2;
              if (!sf_getint("niter",&niter)) niter=100; /*smooth division maximum iterations*/
          } else {
              if (!sf_getfloat("eps", &eps)) eps=SF_EPS; /*padding*/
              if (!sf_getfloat("perc", &perc)) perc=0.1; /*percentage of maximum for padding*/
              perc /= 100;
          }
      }
      if (stable==3) {
        if (!sf_getbool("freq_scal", &freq_scal)) freq_scal=false; /*frequency amplitude spectrum scaling*/
        if (!sf_getfloat("eps", &eps)) eps=SF_EPS; /*padding*/
        if (!sf_getfloat("perc", &perc)) perc=0.1; /*percentage of maximum for padding*/
        perc /= 100;
      }
    }
    /* source/receiver info */
    if (!sf_getint("shtbgn", &shtbgn)) sf_error("Need shot starting location on grid!");
    if (!sf_getint("sht0", &sht0)) sht0=shtbgn; /*actual shot origin on grid*/
    if (!sf_getint("shtend", &shtend)) sf_error("Need shot ending location on grid!");
    if (!sf_getint("shtint", &shtint)) sf_error("Need shot interval on grid!");
    shtnum = (int)((shtend-shtbgn)/shtint) + 1;
    shtnum0 = shtnum;
    if (shtnum%numprocs!=0) {
      shtnum += numprocs-shtnum%numprocs;
      if (verb && cpuid==0) sf_warning("Total shot number is not divisible by total number of nodes! shunum padded from %d to %d.",shtnum0,shtnum); }
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
    /* shot output id */
    if (!sf_getint("shtid", &shtid)) shtid=0;

    /*Set I/O file*/
    if (adj) { /* migration */
      Frcd = sf_input("--input"); /*record from elsewhere*/
      Fsrc  = sf_input("src");   /*source wavelet*/      
      Fimg  = sf_output("--output");
    } else { /* modeling */
      Fimg = sf_input("--input");
      Frcd = sf_output("--output");
      Fsrc  = sf_input("src");   /*source wavelet*/      
    }
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
    /* wavefield axis */
    wfnt = (int)(nt-1)/snpint+1;
    wfdt = dt*snpint;

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
    if (adj) {
      sf_histint(Frcd,"n1", &tmpint);
      if (tmpint != nt ) sf_error("Error parameter n1 in record!");
      sf_histint(Frcd,"n2", &tmpint);
      if (tmpint != gpl ) sf_error("Error parameter n2 in record!");
      sf_histint(Frcd,"n3", &tmpint);
      if (tmpint != shtnum0 ) sf_error("Error parameter n3 in record!");
    }

    /*allocate memory*/
    ww=sf_complexalloc(nt);
    rr=sf_floatalloc(nzx);
    ltf = sf_complexalloc2(nzx,m2);
    rtf = sf_complexalloc2(m2,nk);
    ltb = sf_complexalloc2(nzx,m2b);
    rtb = sf_complexalloc2(m2b,nk);
    geop = (geopar) sf_alloc(1, sizeof(*geop));
    imgsum = sf_complexalloc2(nz, nx);
    if (cpuid==0) {
      record = sf_complexalloc3(nt, gpl, shtnum);
    } else {
      record = NULL;
    }
    if (adj && stable==1) {
      img_visc = sf_complexalloc(nz*nx);
      img_disp = sf_complexalloc(nz*nx);
      ratio  = sf_floatalloc(nz*nx);
      if (sdiv) {
        img_visc_f = sf_floatalloc(nz*nx);
        img_disp_f = sf_floatalloc(nz*nx);
      } else {
        img_visc_f = NULL;
        img_disp_f = NULL;
      }
    } else {
      img_visc = NULL;
      img_disp = NULL;
      ratio  = NULL;
      img_visc_f = NULL;
      img_disp_f = NULL;
    }

    /*read from files*/
    sf_complexread(ww,nt,Fsrc);
    sf_complexread(ltf[0],nzx*m2,left);
    sf_complexread(rtf[0],m2*nk,right);
    sf_complexread(ltb[0],nzx*m2b,leftb);
    sf_complexread(rtb[0],m2b*nk,rightb);
    if(adj) {
      /*read data*/
      if (cpuid==0) {
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
      }
      /*initialize image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
      for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	  imgsum[ix][iz] = sf_cmplx(0.,0.);
    } else {
      /*read image*/
      sf_complexread(imgsum[0],nx*nz,Fimg);
      /*initialize data*/
      if (cpuid==0) {
#ifdef _OPENMP
#pragma omp parallel for private(is,ix,it)
#endif
      for (is=0; is<shtnum; is++)
	for (ix=0; ix<gpl; ix++)
	  for (it=0; it<nt; it++)
	    record[is][ix][it] = sf_cmplx(0.,0.);
      }
    }

    /* output RSF files */
    if (cpuid==0) {
      if (adj) { /* migration */
	sf_setn(ax, nx);
	sf_setn(az, nz);
	/*write image*/
	sf_oaxa(Fimg, az, 1);
	sf_oaxa(Fimg, ax, 2);
	sf_settype(Fimg,SF_COMPLEX);
	if (stable==1 && NULL!=sf_getstring("norm")) {
	  Fnorm = sf_output("norm");
	  sf_oaxa(Fnorm, az, 1);
	  sf_oaxa(Fnorm, ax, 2);
	  sf_settype(Fnorm,SF_FLOAT);
	} else Fnorm = NULL;
      } else { /* modeling */
	sf_setn(ax, gpl);
	as = sf_iaxa(Fvel, 2);
	sf_setn(as,shtnum0);
	sf_setd(as,shtint*dx);
	sf_seto(as,shtbgn*dx+ox);
	sf_oaxa(Frcd, at, 1);
	sf_oaxa(Frcd, ax, 2);
	sf_oaxa(Frcd, as ,3);
	sf_settype(Frcd,SF_COMPLEX);
      }
    }

    if (cpuid == shtid%numprocs) {
      if (NULL!=sf_getstring("tmpwf")) {
	Ftmpwf  = sf_output("tmpwf"); /*wavefield snap*/
	sf_setn(ax, nx);
	sf_setn(az, nz);
      	sf_setn(at, wfnt);
	sf_setd(at, wfdt);
	sf_oaxa(Ftmpwf, az, 1);
	sf_oaxa(Ftmpwf, ax, 2);
	sf_oaxa(Ftmpwf, at, 3);
	sf_settype(Ftmpwf,SF_COMPLEX);
      } else Ftmpwf = NULL;

      if (NULL!=sf_getstring("tmpwfb")) {
	Ftmpwfb  = sf_output("tmpwfb"); /*wavefield snap*/
	sf_setn(ax, nx);
	sf_setn(az, nz);
      	sf_setn(at, wfnt);
	sf_setd(at, wfdt);
	sf_oaxa(Ftmpwfb, az, 1);
	sf_oaxa(Ftmpwfb, ax, 2);
	sf_oaxa(Ftmpwfb, at, 3);
	sf_settype(Ftmpwfb,SF_COMPLEX);
      } else Ftmpwfb = NULL;
    } else {
      Ftmpwf = NULL;
      Ftmpwfb = NULL;
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
    geop->adj   = adj;
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
    
    if (adj) {
      switch (stable) {
      case 0:
	/*conventional imaging condition*/
        if (NULL == Ftmpwf && NULL == Ftmpwfb)
          psrtm(record, imgsum, geop);
        else
          psrtm_mov(record, imgsum, geop, Ftmpwf, Ftmpwfb, shtid);
	break;

      case 1:
	/*stable imaging condition for Q-compensation*/
	geop->mode = 0;	/*preceding the user specification*/
	psrtm(record, imgsum, geop);
	/*save the viscoacoustic image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    img_visc[iz+ix*nz] = imgsum[ix][iz];
          }
        }

	geop->mode = 1;
	psrtm(record, imgsum, geop);
	/*save the dispersion-only image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    img_disp[iz+ix*nz] = imgsum[ix][iz];
          }
        }

        /* stable division */
        if (sdiv) {
	  /*smooth division*/
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
          for (iz=0; iz<nz*nx; iz++) {
	      img_visc_f[iz] = cabsf(img_visc[iz]);
	      img_disp_f[iz] = cabsf(img_disp[iz]);
          }
	  dims[0] = nz; dims[1] = nx;
	  sf_divn_init(2, nx*nz, dims, rect, niter, verb);
	  sf_divne (img_disp_f, img_visc_f, ratio, reg);
        } else {
          vmax = find_cmax(nz*nx, img_visc);
          veps = (vmax==0) ? eps : vmax*vmax*perc;
          if (verb) sf_warning("Space domain: vmax=%f, veps=%f",vmax,veps);
          stable_cdiv_f(nz*nx, veps, img_disp, img_visc, ratio);
        }
	  
	/*apply stable imaging condition*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
            imgsum[ix][iz] = img_disp[iz+ix*nz]*ratio[iz+ix*nz];
#else
            imgsum[ix][iz] = sf_crmul(img_disp[iz+ix*nz],ratio[iz+ix*nz]);
#endif
	  }
	}
	break;
	
      case 2: 
	/*shot-by-shot q-rtm with normalization*/
	psqrtm_sbs(record, imgsum, geop, sdiv);
	break;

      case 3:
	/*snapshot-by-snapshot compensation*/
	psqrtm_com(record, imgsum, geop, freq_scal, Ftmpwf, Ftmpwfb, shtid);
	break;

      case 4:
	/*deconvolution imaging condition*/
	psqrtm_dec(record, imgsum, geop);
	break;

      default:
	sf_error("unknown stable %d",stable);
      }

    } else {
      /*born modelling*/
      psrtm(record, imgsum, geop);
    }
    
    if (cpuid==0) {
      if (adj) {
	sf_complexwrite(imgsum[0], nx*nz, Fimg);
	if (NULL != Fnorm) sf_floatwrite(ratio, nx*nz, Fnorm);
      } else {
	sf_complexwrite(record[0][0], shtnum0*gpl*nt, Frcd);
      }
    }

    /*free memory*/
    free(geop);
    free(ww); free(rr);
    free(*ltf); free(ltf);
    free(*rtf); free(rtf);
    free(*ltb); free(ltb);
    free(*rtb); free(rtb);
    free(*imgsum); free(imgsum);
    if (adj && stable==1) { 
      free(img_visc); free(img_disp); free(ratio);
    }
    if (cpuid==0) {
      free(**record); free(*record); free(record);
    }

    MPI_Finalize();
    exit(0);
}

int psrtm(sf_complex*** record, sf_complex** imgsum, geopar geop)
/*< low-rank one-step pre-stack RTM linear operator >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool adj; /* migration(adjoint) flag */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
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
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    m2 = geop->m2; m2b = geop->m2b; pad1 = geop->pad1;
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

      if(adj) {
	if (cpuid==0) sendbuf = record[is*numprocs][0];
	else sendbuf = NULL;
	recvbuf = tmprec[0];
	MPI_Scatter(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // tmprec[ix][it] = record[is][ix][it];
      }
      
      if (shtcur<shtnum0) {
	lrosback2(img, wvfld, sill, tmprec, geop);
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
	if (cpuid==0) recvbuf = record[is*numprocs][0];
	else recvbuf = NULL;
	sendbuf = tmprec[0];
	MPI_Gather(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // record[is][ix][it] = tmprec[ix][it];
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

int psqrtm_sbs(sf_complex*** record, sf_complex** imgsum, geopar geop, bool sdiv)
/*< low-rank one-step pre-stack RTM linear operator >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool adj; /* migration(adjoint) flag -> not used in this case */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    /*pointers*/
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    sf_complex ***wvfld;
    sf_complex **tmprec, **rec, **img;
    float **sill;
    
    /*normolization*/
    sf_complex *img_visc, *img_disp;
    float *ratio;
    float *img_visc_f, *img_disp_f;

    /*misc*/
    int ix, iz, is;
    int wfnt;
    float wfdt;
    float vmax,veps;

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
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    m2 = geop->m2; m2b = geop->m2b; pad1 = geop->pad1;
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
    rec = sf_complexalloc2(nt, gpl);
    img_visc = sf_complexalloc(nz*nx);
    img_disp = sf_complexalloc(nz*nx);
    ratio    = sf_floatalloc(nz*nx);
    if (sdiv) {
        img_visc_f = sf_floatalloc(nz*nx);
        img_disp_f = sf_floatalloc(nz*nx);
    } else {
        img_visc_f = NULL;
        img_disp_f = NULL;
    }

#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
      for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	  imgsum[ix][iz] = sf_cmplx(0.,0.);

    /* start the work */
    tstart = clock();

    for (is=0; is*numprocs<shtnum; is++){

      shtcur = is*numprocs+cpuid; // current shot index

      /*scatter the data*/
      if (cpuid==0) sendbuf = record[is*numprocs][0];
      else sendbuf = NULL;
      recvbuf = rec[0];
      MPI_Scatter(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // rec[ix][it] = record[is][ix][it];

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

	geop->mode = 0; /*preceding the user specification -> viscoacoustic*/
	lrosfor2(wvfld, sill, tmprec, geop);
	lrosback2(img, wvfld, sill, rec, geop);

	/*save the viscoacoustic image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    img_visc[iz+ix*nz] = img[ix][iz];
          }
        }

	geop->mode = 1; /*preceding the user specification -> dispersion-only*/
	lrosfor2(wvfld, sill, tmprec, geop);
	lrosback2(img, wvfld, sill, rec, geop);

	/*save the dispersion-only image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    img_disp[iz+ix*nz] = img[ix][iz];
          }
        }

        /* stable division */
        if (sdiv) {
          /*smooth division*/
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
          for (iz=0; iz<nz*nx; iz++) {
            img_visc_f[iz] = cabsf(img_visc[iz]);
            img_disp_f[iz] = cabsf(img_disp[iz]);
          }
          dims[0] = nz; dims[1] = nx;
          sf_divn_init(2, nx*nz, dims, rect, niter, verb);
          sf_divne (img_disp_f, img_visc_f, ratio, reg);
        } else {
          vmax = find_cmax(nz*nx, img_visc);
          veps = (vmax==0) ? eps : vmax*vmax*perc;
          if (verb) sf_warning("Space domain: vmax=%f, veps=%f",vmax,veps);
          stable_cdiv_f(nz*nx, veps, img_disp, img_visc, ratio);
        }

	/*apply stable imaging condition*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
            img[ix][iz] = img_disp[iz+ix*nz]*ratio[iz+ix*nz];
#else
            img[ix][iz] = sf_crmul(img_disp[iz+ix*nz],ratio[iz+ix*nz]);
#endif
	  }
	}

	/* image reduction */
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

      } /*if (shtcur<shtnum0)*/

    } /*shot iteration*/

    MPI_Barrier(MPI_COMM_WORLD);

    /*image reduction*/
#if MPI_VERSION >= 2
    sendbuf = (sf_complex *) MPI_IN_PLACE;
#else /* will fail */
    sendbuf = NULL;
#endif 
    recvbuf = imgsum[0];
    MPI_Allreduce(sendbuf, recvbuf, nx*nz, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    if (verb)
      sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

    free(**wvfld); free(*wvfld); free(wvfld);
    free(*tmprec); free(tmprec);
    free(*rec); free(rec);
    free(*img); free(img);
    free(img_visc); free(img_disp); free(ratio);
    if (illum) { free(*sill); free(sill);}

    return 0;
}

int psqrtm_com(sf_complex*** record, sf_complex** imgsum, geopar geop, bool freq_scal, sf_file Ftmpwf, sf_file Ftmpwfb, int shtid)
/*< low-rank one-step pre-stack RTM linear operator >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool adj; /* migration(adjoint) flag -> not used in this case */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    /*pointers*/
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    sf_complex ***wvfld, ***wvfld_b;
    sf_complex **tmprec, **rec, **img;
    float **sill;
    
    /*normolization*/
    sf_complex *wvfld_visc, *wvfld_disp;
    float *ratio;

    /*frequency domain compensation*/
    kiss_fft_cpx **ctrace3=NULL;
    sf_complex **trace3=NULL, *cc=NULL, *tmp_inp=NULL;
    float *wvfld_visc_f=NULL, *wvfld_disp_f=NULL, *ratio_f=NULL, *tmp_out=NULL;
    kiss_fft_cfg *cfg3=NULL, *icfg3=NULL;
    int n3,n2,n1,i3,i2,i1,nth;
    float wt;

    /*misc*/
    int ix, iz, it, is, ith;
    int wfnt;
    float wfdt;
    float vmax,veps;

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
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    m2 = geop->m2; m2b = geop->m2b; pad1 = geop->pad1;
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
    wvfld_b = sf_complexalloc3(nz, nx, wfnt);
    if (illum) sill = sf_floatalloc2(nz, nx);
    else sill = NULL;
    img = sf_complexalloc2(nz, nx);
    rec = sf_complexalloc2(nt, gpl);
    wvfld_visc = sf_complexalloc(nz*nx*wfnt);
    wvfld_disp = sf_complexalloc(nz*nx*wfnt);
    ratio      = sf_floatalloc(nz*nx*wfnt);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
      for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	  imgsum[ix][iz] = sf_cmplx(0.,0.);

    if (freq_scal) {
        n1 = nz;
        n2 = nx;
        n3 = kiss_fft_next_fast_size(nt);
        wvfld_visc_f = sf_floatalloc(n1*n2*n3);
        wvfld_disp_f = sf_floatalloc(n1*n2*n3);
        ratio_f      = sf_floatalloc(n1*n2*n3);
        cc = sf_complexalloc(n1*n2*n3);
#ifdef _OPENMP
#pragma omp parallel
        {nth = omp_get_num_threads();}
#endif
        cfg3  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
        icfg3 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
        for (ith=0; ith<nth; ith++) {
            cfg3[ith] = kiss_fft_alloc(n3,0,NULL,NULL);
            icfg3[ith]= kiss_fft_alloc(n3,1,NULL,NULL);
        }
        trace3 = sf_complexalloc2(n3,nth);
        ctrace3= (kiss_fft_cpx **) trace3;
        wt =  1.0/n3;
    }

    /* start the work */
    tstart = clock();

    for (is=0; is*numprocs<shtnum; is++){

      shtcur = is*numprocs+cpuid; // current shot index

      /*scatter the data*/
      if (cpuid==0) sendbuf = record[is*numprocs][0];
      else sendbuf = NULL;
      recvbuf = rec[0];
      MPI_Scatter(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // rec[ix][it] = record[is][ix][it];

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

	/*source wavefield*/
	geop->mode = 0; /*preceding the user specification -> viscoacoustic*/
	lrosfor2q(wvfld, geop);
	/*save the viscoacoustic source wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
	for (it=0; it<wfnt; it++) {
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++){
	      wvfld_visc[iz+ix*nz+it*nx*nz] = wvfld[it][ix][iz];
            }
          }
        }

	geop->mode = 1; /*preceding the user specification -> dispersion-only*/
	if (illum) { /*illumination by source wavefield*/
	  lrosfor2(wvfld, sill, tmprec, geop);
	} else {
	  lrosfor2q(wvfld, geop);
	}
	/*save the dispersion-only source wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
	for (it=0; it<wfnt; it++) {
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
	      wvfld_disp[iz+ix*nz+it*nx*nz] = wvfld[it][ix][iz];
            }
          }
        }

        /* stable division */
        if (freq_scal) {

            /***************************************/
            /* visc -> visc_f */
            tmp_inp = wvfld_visc;
            tmp_out = wvfld_visc_f;

            /********* Forward FFT **********/

#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (i3<nt) {
                            cc[(i3*n2+i2)*n1+i1] = tmp_inp[(i3*nx+i2)*nz+i1];
                        } else {
                            cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                        }
                    }
                }
            }

            /* FFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(cfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
#ifdef SF_HAS_COMPLEX_H
                        tmp_out[(i3*n2+i2)*n1+i1] = cabsf(trace3[ith][i3]);
#else
                        tmp_out[(i3*n2+i2)*n1+i1] = sf_cabsf(trace3[ith][i3]);
#endif
                    }
                }
            }

            /***************************************/
            /* disp -> disp_f */
            tmp_inp = wvfld_disp;
            tmp_out = wvfld_disp_f;

            /********* Forward FFT **********/

#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (i3<nt) {
                            cc[(i3*n2+i2)*n1+i1] = tmp_inp[(i3*nx+i2)*nz+i1];
                        } else {
                            cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                        }
                    }
                }
            }

            /* FFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(cfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
                        cc[(i3*n2+i2)*n1+i1] = trace3[ith][i3];
#ifdef SF_HAS_COMPLEX_H
                        tmp_out[(i3*n2+i2)*n1+i1] = cabsf(trace3[ith][i3]);
#else
                        tmp_out[(i3*n2+i2)*n1+i1] = sf_cabsf(trace3[ith][i3]);
#endif
                    }
                }
            }

            /* Stable division of amplitude spectrum */
            //for (it=0; it<n3; it++) {
            //    vmax = find_max(nz*nx, wvfld_visc_f+it*nz*nx);
            //    veps = (vmax==0) ? eps : vmax*vmax*perc;
            //    if (verb) sf_warning("Frequency domain: vmax=%f, veps=%f",vmax,veps);
            //    stable_div(nz*nx, veps, wvfld_disp_f+it*nz*nx, wvfld_visc_f+it*nz*nx, ratio_f+it*nz*nx);
            //}
            vmax = find_max(n1*n2*n3, wvfld_visc_f);
            veps = (vmax==0) ? eps : vmax*vmax*perc;
            if (verb) sf_warning("Frequency domain: vmax=%f, veps=%f",vmax,veps);
            stable_div(n1*n2*n3, veps, wvfld_disp_f, wvfld_visc_f, ratio_f);

            /* Scaling the wavefield */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        cc[(i3*n2+i2)*n1+i1] = cc[(i3*n2+i2)*n1+i1]*ratio_f[(i3*n2+i2)*n1+i1];
#else
                        cc[(i3*n2+i2)*n1+i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],ratio_f[(i3*n2+i2)*n1+i1]);
#endif
                    }
                }
            }

            /********* Inverse FFT **********/

            /* IFFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(icfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
                        cc[(i3*n2+i2)*n1+i1] = trace3[ith][i3];
                    }
                }
            }

            /* FFT normalization */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<nt; i3++) {
                for (i2=0; i2<nx; i2++) {
                    for (i1=0; i1<nz; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        wvfld[i3][i2][i1] = wt*cc[(i3*n2+i2)*n1+i1];
#else
                        wvfld[i3][i2][i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],wt);
#endif
                    }
                }
            }

        } else {
            for (it=0; it<wfnt; it++) {
                vmax = find_cmax(nz*nx, wvfld_visc+it*nz*nx);
                veps = (vmax==0) ? eps : vmax*vmax*perc;
                if (verb) sf_warning("Space domain: vmax=%f, veps=%f",vmax,veps);
                stable_cdiv_f(nz*nx, veps, wvfld_disp+it*nz*nx, wvfld_visc+it*nz*nx, ratio+it*nz*nx);
            }

            /*scaling the source wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
            for (it=0; it<wfnt; it++) {
                for (ix=0; ix<nx; ix++) {
                    for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
                        wvfld[it][ix][iz] = wvfld_disp[iz+ix*nz+it*nx*nz]*ratio[iz+ix*nz+it*nx*nz];
#else
                        wvfld[it][ix][iz] = sf_crmul(wvfld_disp[iz+ix*nz+it*nx*nz],ratio[iz+ix*nz+it*nx*nz]);
#endif
                    }
                }
            }

        }

	/*receiver wavefield*/
	geop->mode = 0; /*preceding the user specification -> viscoacoustic*/
	lrosback2q(wvfld_b, rec, geop);
	/*save the viscoacoustic receiver wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
	for (it=0; it<wfnt; it++) {
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++){
	      wvfld_visc[iz+ix*nz+it*nx*nz] = wvfld_b[it][ix][iz];
            }
          }
        }

	geop->mode = 1; /*preceding the user specification -> dispersion-only*/
	lrosback2q(wvfld_b, rec, geop);
	/*save the dispersion-only receiver wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
	for (it=0; it<wfnt; it++) {
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
	      wvfld_disp[iz+ix*nz+it*nx*nz] = wvfld_b[it][ix][iz];
            }
          }
        }

        /* stable division */
        if (freq_scal) {

            /***************************************/
            /* visc -> visc_f */
            tmp_inp = wvfld_visc;
            tmp_out = wvfld_visc_f;

            /********* Forward FFT **********/

#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (i3<nt) {
                            cc[(i3*n2+i2)*n1+i1] = tmp_inp[(i3*nx+i2)*nz+i1];
                        } else {
                            cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                        }
                    }
                }
            }

            /* FFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(cfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
#ifdef SF_HAS_COMPLEX_H
                        tmp_out[(i3*n2+i2)*n1+i1] = cabsf(trace3[ith][i3]);
#else
                        tmp_out[(i3*n2+i2)*n1+i1] = sf_cabsf(trace3[ith][i3]);
#endif
                    }
                }
            }

            /***************************************/
            /* disp -> disp_f */
            tmp_inp = wvfld_disp;
            tmp_out = wvfld_disp_f;

            /********* Forward FFT **********/

#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (i3<nt) {
                            cc[(i3*n2+i2)*n1+i1] = tmp_inp[(i3*nx+i2)*nz+i1];
                        } else {
                            cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                        }
                    }
                }
            }

            /* FFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(cfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
                        cc[(i3*n2+i2)*n1+i1] = trace3[ith][i3];
#ifdef SF_HAS_COMPLEX_H
                        tmp_out[(i3*n2+i2)*n1+i1] = cabsf(trace3[ith][i3]);
#else
                        tmp_out[(i3*n2+i2)*n1+i1] = sf_cabsf(trace3[ith][i3]);
#endif
                    }
                }
            }

            /* Stable division of amplitude spectrum */
            //for (it=0; it<n3; it++) {
            //    vmax = find_max(nz*nx, wvfld_visc_f+it*nz*nx);
            //    veps = (vmax==0) ? eps : vmax*vmax*perc;
            //    if (verb) sf_warning("Frequency domain: vmax=%f, veps=%f",vmax,veps);
            //    stable_div(nz*nx, veps, wvfld_disp_f+it*nz*nx, wvfld_visc_f+it*nz*nx, ratio_f+it*nz*nx);
            //}
            vmax = find_max(n1*n2*n3, wvfld_visc_f);
            veps = (vmax==0) ? eps : vmax*vmax*perc;
            if (verb) sf_warning("Frequency domain: vmax=%f, veps=%f",vmax,veps);
            stable_div(n1*n2*n3, veps, wvfld_disp_f, wvfld_visc_f, ratio_f);

            /* Scaling the wavefield */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        cc[(i3*n2+i2)*n1+i1] = cc[(i3*n2+i2)*n1+i1]*ratio_f[(i3*n2+i2)*n1+i1];
#else
                        cc[(i3*n2+i2)*n1+i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],ratio_f[(i3*n2+i2)*n1+i1]);
#endif
                    }
                }
            }

            /********* Inverse FFT **********/

            /* IFFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                for (i1=0; i1<n1; i1++) {
                    kiss_fft_stride(icfg3[ith],(kiss_fft_cpx *) (cc+i2*n1+i1),ctrace3[ith],n2*n1);
                    for (i3=0; i3<n3; i3++) {
                        cc[(i3*n2+i2)*n1+i1] = trace3[ith][i3];
                    }
                }
            }

            /* FFT normalization */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<nt; i3++) {
                for (i2=0; i2<nx; i2++) {
                    for (i1=0; i1<nz; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        wvfld_b[i3][i2][i1] = wt*cc[(i3*n2+i2)*n1+i1];
#else
                        wvfld_b[i3][i2][i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],wt);
#endif
                    }
                }
            }

        } else {
            for (it=0; it<wfnt; it++) {
                vmax = find_cmax(nz*nx, wvfld_visc+it*nz*nx);
                veps = (vmax==0) ? eps : vmax*vmax*perc;
                if (verb) sf_warning("Space domain: vmax=%f, veps=%f",vmax,veps);
                stable_cdiv_f(nz*nx, veps, wvfld_disp+it*nz*nx, wvfld_visc+it*nz*nx, ratio+it*nz*nx);
            }

            /*scaling the receiver wavefield*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
            for (it=0; it<wfnt; it++) {
                for (ix=0; ix<nx; ix++) {
                    for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
                        wvfld_b[it][ix][iz] = wvfld_disp[iz+ix*nz+it*nx*nz]*ratio[iz+ix*nz+it*nx*nz];
#else
                        wvfld_b[it][ix][iz] = sf_crmul(wvfld_disp[iz+ix*nz+it*nx*nz],ratio[iz+ix*nz+it*nx*nz]);
#endif
                    }
                }
            }
        }

	/*cross-correlation imaging condition*/
	ccrimg(img, wvfld, wvfld_b, sill, geop);

	/* image reduction */
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

	if (NULL!=Ftmpwf && shtcur==shtid) {
          sf_warning("Writing source wavefield on cpu %d",cpuid);
	  sf_complexwrite(wvfld[0][0], wfnt*nx*nz, Ftmpwf);
        }
	if (NULL!=Ftmpwfb && shtcur==shtid) {
          sf_warning("Writing receiver wavefield on cpu %d",cpuid);
	  sf_complexwrite(wvfld_b[0][0], wfnt*nx*nz, Ftmpwfb);
        }

      } /*if (shtcur<shtnum0)*/

    } /*shot iteration*/

    MPI_Barrier(MPI_COMM_WORLD);

    /*image reduction*/
#if MPI_VERSION >= 2
    sendbuf = (sf_complex *) MPI_IN_PLACE;
#else /* will fail */
    sendbuf = NULL;
#endif 
    recvbuf = imgsum[0];
    MPI_Allreduce(sendbuf, recvbuf, nx*nz, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    if (verb)
      sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

    free(**wvfld); free(*wvfld); free(wvfld);
    free(**wvfld_b); free(*wvfld_b); free(wvfld_b);
    free(*tmprec); free(tmprec);
    free(*rec); free(rec);
    free(*img); free(img);
    free(wvfld_visc); free(wvfld_disp); free(ratio);
    if (illum) { free(*sill); free(sill);}
    if (freq_scal) {
        for (ith=0; ith<nth; ith++) {
            if (NULL!=cfg3)  { free(cfg3[ith]);   cfg3[ith]=NULL; }
            if (NULL!=icfg3) { free(icfg3[ith]); icfg3[ith]=NULL; }
        }
        if (NULL!=cc) { free(cc); cc=NULL; }
        if (NULL!=trace3) { free(*trace3); free(trace3); trace3=NULL; }
    }

    return 0;
}

int psqrtm_dec(sf_complex*** record, sf_complex** imgsum, geopar geop)
/*< low-rank one-step pre-stack RTM linear operator >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool adj; /* migration(adjoint) flag -> not used in this case */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    /*pointers*/
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    sf_complex ***wvfld;
    sf_complex **tmprec, **rec, **img;
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
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    m2 = geop->m2; m2b = geop->m2b; pad1 = geop->pad1;
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
    sill = sf_floatalloc2(nz, nx);
    img = sf_complexalloc2(nz, nx);
    rec = sf_complexalloc2(nt, gpl);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
      for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	  imgsum[ix][iz] = sf_cmplx(0.,0.);

    /* start the work */
    tstart = clock();

    for (is=0; is*numprocs<shtnum; is++){

      shtcur = is*numprocs+cpuid; // current shot index

      /*scatter the data*/
      if (cpuid==0) sendbuf = record[is*numprocs][0];
      else sendbuf = NULL;
      recvbuf = rec[0];
      MPI_Scatter(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // rec[ix][it] = record[is][ix][it];

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
	  sf_warning("Using receiver illumination for deconvolution imaging condition!");
	  sf_warning("============================");
	}
	
	/*generate reflectivity map*/
	reflgen(nzb, nxb, spz+top, spx+lft, rectz, rectx, repeat, rr);

	geop->mode = 1;/*dispersion-only*/
	rcvill2(sill, rec, geop);

	geop->mode = 0; /*viscoacoustic*/
	geop->illum= false;
	lrosfor2(wvfld, sill, tmprec, geop);
	geop->illum= true;
	lrosback2(img, wvfld, sill, rec, geop);

	/* image reduction */
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

      } /*if (shtcur<shtnum0)*/

    } /*shot iteration*/

    MPI_Barrier(MPI_COMM_WORLD);

    /*image reduction*/
#if MPI_VERSION >= 2
    sendbuf = (sf_complex *) MPI_IN_PLACE;
#else /* will fail */
    sendbuf = NULL;
#endif 
    recvbuf = imgsum[0];
    MPI_Allreduce(sendbuf, recvbuf, nx*nz, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    if (verb)
      sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

    free(**wvfld); free(*wvfld); free(wvfld);
    free(*tmprec); free(tmprec);
    free(*rec); free(rec);
    free(*img); free(img);
    free(*sill); free(sill);

    return 0;
}

int psrtm_mov(sf_complex*** record, sf_complex** imgsum, geopar geop, sf_file Ftmpwf, sf_file Ftmpwfb, int shtid)
/*< low-rank one-step pre-stack RTM linear operator that outputs movie >*/
{
    /*geopar variables*/
    int nx, nz;
    int nxb, nzb;
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl; /*source/geophone location*/
    int snpint;
    int top, bot, lft, rht; /*abc boundary*/
    int nt;
    float dt;
    float trunc; 
    bool adj; /* migration(adjoint) flag -> not used in this case */
    bool verb; /* verbosity flag */
    bool illum; /* source illumination flag*/
    int m2, m2b, pad1;
    /*pointers*/
    float *rr;
    /*extras*/
    bool roll; /* survey strategy */
    int rectz,rectx,repeat; /*refl smoothing parameters*/
    int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
    /*mpi*/
    int cpuid, numprocs;

    sf_complex ***wvfld, ***wvfld_b;
    sf_complex **tmprec, **rec, **img;
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
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    m2 = geop->m2; m2b = geop->m2b; pad1 = geop->pad1;
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
    wvfld_b = sf_complexalloc3(nz, nx, wfnt);
    if (illum) sill = sf_floatalloc2(nz, nx);
    else sill = NULL;
    img = sf_complexalloc2(nz, nx);
    rec = sf_complexalloc2(nt, gpl);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
      for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	  imgsum[ix][iz] = sf_cmplx(0.,0.);

    /* start the work */
    tstart = clock();

    for (is=0; is*numprocs<shtnum; is++){

      shtcur = is*numprocs+cpuid; // current shot index

      /*scatter the data*/
      if (cpuid==0) sendbuf = record[is*numprocs][0];
      else sendbuf = NULL;
      recvbuf = rec[0];
      MPI_Scatter(sendbuf, gpl*nt, MPI_COMPLEX, recvbuf, gpl*nt, MPI_COMPLEX, 0, MPI_COMM_WORLD); // rec[ix][it] = record[is][ix][it];

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

	/*source wavefield*/
        lrosfor2(wvfld, sill, tmprec, geop);

	/*receiver wavefield*/
	lrosback2q(wvfld_b, rec, geop);

	/*cross-correlation imaging condition*/
	ccrimg(img, wvfld, wvfld_b, sill, geop);

	/* image reduction */
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

	if (NULL!=Ftmpwf && shtcur==shtid)
	  sf_complexwrite(wvfld[0][0], wfnt*nx*nz, Ftmpwf);
	if (NULL!=Ftmpwfb && shtcur==shtid)
	  sf_complexwrite(wvfld_b[0][0], wfnt*nx*nz, Ftmpwfb);

      } /*if (shtcur<shtnum0)*/

    } /*shot iteration*/

    MPI_Barrier(MPI_COMM_WORLD);

    /*image reduction*/
#if MPI_VERSION >= 2
    sendbuf = (sf_complex *) MPI_IN_PLACE;
#else /* will fail */
    sendbuf = NULL;
#endif 
    recvbuf = imgsum[0];
    MPI_Allreduce(sendbuf, recvbuf, nx*nz, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    if (verb)
      sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

    free(**wvfld); free(*wvfld); free(wvfld);
    free(**wvfld_b); free(*wvfld_b); free(wvfld_b);
    free(*tmprec); free(tmprec);
    free(*rec); free(rec);
    free(*img); free(img);
    if (illum) { free(*sill); free(sill);}

    return 0;
}
