/* 2-D two-components elastic wavefield modeling operators with lowrank approximation. */
/*
  Copyright (C) 2014 UT Austin
     
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
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

/*automatically generated headers*/
#include "muting.h"
#include "absorb.h"
#include "fft2w.h"
#include "puthead.h"
#include "ricker1.h"

/* wave propagation utils*/
typedef struct Geopar {
  /*geometry*/
  int   nx, nz;
  int   nx1, nz1;
  float dx, dz;
  int   gpz, gpx, gpl;
  /*source parameters*/
  int nt;
  float dt;
  int src;
  int isx, isz;
  float f0, t0, A;
  /*absorbing boundary*/
  bool abc;
  int nbt, nbb, nbl, nbr; /* boundary width */
  float ct,cb,cl,cr; /* decay parameter */
  /*propagator*/
  int m2p, m2s;
  float **plt, **prt, **slt, **srt;
  /*wavefield*/
  int snpint, wfnt;
  /*verbosity*/
  bool verb;
  /*fft parameters*/
  bool cmplx;
  int pad1;
  int nk;
} * geopar; /*geometry parameters*/

typedef struct Mpipar {
    int cpuid;
    int numprocs;
} * mpipar; /*geometry parameters*/

int lrewefor2(float **rcdz, float** rcdx, float **upz, float **upx, float **usz, float **usx, geopar geop);
int lrewebac2(float *imgpp, float *imgps, float *imgsp, float *imgss, float **rcdz, float** rcdx, float **upz, float **upx, float **usz, float **usx, geopar geop);


/**************************************************************************/
/* MAIN FUNC */
int main(int argc, char* argv[])
{
  /* modeling parameters (passed to engine) */
  int nx,nz,nz1,nx1;     /* model dimension parameters */
  float dx,dz;           /* spatial sampling */
  int gpz;               /* geophone depth */

  int nt, src, isx, isz; /* model/data dimension and source location */
  float dt, f0, t0, A;   /* time and source related */

  bool abc;              /* abc flag */
  int nbt=0,nbb=0,nbl=0,nbr=0;   /* abc width */
  float ct,cb,cl,cr;     /* decay parameter */

  int m2p,m2s;           /* lowrank matrices dimension */
  float **plt, **prt, **slt, **srt; /* lowrank matrices */

  int snpint, wfnt; /* snapshot interval and total number */

  bool verb;  /* verbose flag */
  bool cmplx; /* complex fft flag */
  int pad1;   /* fft padding on the first axis */

  /* handling parameters (used within main function) */
  bool mig;

  geopar geop;

  float wfdt;
  int ix,it,is;
  int nzx,nzx1,nzx2,nz2,nx2,nk,n2; /*fft related*/
  int nth;

  float **rcdz, **rcdx;
  float *imgpp, *imgps, *imgsp, *imgss;
  float **wavepz, **wavepx, **wavesz, **wavesx; /* wavefield array */
  float x0, z0; /* starting location */

  sf_file Fvp0, Fpl, Fpr, Fsl, Fsr; /* input files, lowrank matrices */
  sf_file Fo1, Fo2; /* output files */
  sf_file Fpz, Fpx, Fsz, Fsx;
  sf_file Fimgpp, Fimgps, Fimgsp, Fimgss;
  sf_axis az, ax; /* axes */

  /* shot geometry */
  int shtbgn, shtend, shtint;
  int shtnum, shtnum0, shtcur;

  /* muting */
  bool mute;
  int wd;
  float vref;

  /* timing */
  clock_t tstart,tend;
  double duration;

  /* MPI */
  int rank, nodes;
  float *sendbuf, *recvbuf;

  /* I/O arrays */
  float ***recordz, ***recordx;
  /*float *imgpp_sum, *imgps_sum, *imgsp_sum, *imgss_sum;*/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);

  sf_init(argc,argv);

  if(rank==0) sf_warning("nodes=%d",nodes);

  /***************************************************************/
  /* setup Input files */

  Fvp0 = sf_input ("input");  /* vp0 using standard input, just to get model dimension */
  Fpl = sf_input ("pleft");   /* left matrix for P */
  Fpr = sf_input ("pright");  /* right matrix for P */
  Fsl = sf_input ("sleft");   /* left matrix for Sv */
  Fsr = sf_input ("sright");  /* right matrix for Sv */

  /* Read/Write axes */
  az = sf_iaxa(Fvp0,1); nz = sf_n(az); dz = sf_d(az);
  ax = sf_iaxa(Fvp0,2); nx = sf_n(ax); dx = sf_d(ax);
  x0=sf_o(ax);
  z0=sf_o(az);

  /***************************************************************/
  /* read paramters for the commandline */

  /* migration or modeling */
  if(!sf_getbool("mig",&mig)) mig=false; /* migration flag */
  /* time sampling */
  if (!sf_getint("nt",&nt)) nt=301;
  if (!sf_getfloat("dt",&dt)) dt=0.001;
  if (!sf_getint("snpint",&snpint)) snpint=10;
  /* absorbing boundary condition*/
  if(!sf_getbool("abc",&abc)) abc=false; /* absorbing flag */
  if (abc) {
    if(!sf_getint("nbt",&nbt)) sf_error("Need nbt!");
    if(!sf_getint("nbb",&nbb)) nbb = nbt;
    if(!sf_getint("nbl",&nbl)) nbl = nbt;
    if(!sf_getint("nbr",&nbr)) nbr = nbt;
    if(!sf_getfloat("ct",&ct)) sf_error("Need ct!");
    if(!sf_getfloat("cb",&cb)) cb = ct;
    if(!sf_getfloat("cl",&cl)) cl = ct;
    if(!sf_getfloat("cr",&cr)) cr = ct;
  }
  /* source definition */
/*if (!sf_getint("isx",&isx)) isx=nx/2;*/
  if (!sf_getint("isz",&isz)) sf_error("Need shot depth!");
  if (!sf_getfloat("t0",&t0)) t0=0.04;   /*  wavelet time lag */
  if (!sf_getfloat("f0",&f0)) f0=30.0;   /*  wavelet peak freq */
  if (!sf_getfloat("A",&A)) A=1.0;       /*  wavelet amplitude */
  if (!sf_getint("src",&src)) src=1; /* source mode: 1 - exploding force; 2 - equil-energy force */
  /* shot geometry */
  if (!sf_getint("shtbgn", &shtbgn)) sf_error("Need shot starting location on grid!");
  if (!sf_getint("shtend", &shtend)) sf_error("Need shot ending location on grid!");
  if (!sf_getint("shtint", &shtint)) sf_error("Need shot interval on grid!");
  shtnum = (int)((shtend-shtbgn)/shtint) + 1;
  shtnum0 = shtnum; /* remember the actual number of shots before padding */
  /* receiver definition */
  if (!sf_getint("gpz",&gpz)) gpz=nbt+5; /* geophone depth */
  /* flags and misc */
  if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
  if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
  if (!sf_getbool("verb",&verb)) verb=false; /* padding factor on the first axis */
  if (!sf_getbool("mute",&mute)) mute=false; /* muting first arrival */
  if (mute) {
    if (!sf_getint("wd",&wd)) wd=5; /* muting width */
    if (!sf_getfloat("vref",&vref)) vref=1500; /* water velocity */
  }


  /***************************************************************/
  /* set up additional dimension parameters */

  /*nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);*/
  if (cmplx) {
    nk = nz2 = kiss_fft_next_fast_size(nz*pad1);
  } else {
    nk = kiss_fft_next_fast_size(pad1*(nz+1)/2)+1;
    nz2 = 2*(nk-1);
  }
  nx2 = kiss_fft_next_fast_size(nx);
  nk *= nx2;

  /* nz1 and nx1 is the dimension of interest */
  nz1 = nz-nbt-nbb;
  nx1 = nx-nbl-nbr;

  nzx = nz*nx;
  nzx1 = nz1*nx1;
  nzx2 = nz2*nx2;

  wfnt = (int) (nt-1)/snpint + 1;
  wfdt = dt*snpint;


  /***************************************************************/
  /* padding the shot */

  if (shtnum%nodes!=0)
    shtnum += nodes-shtnum%nodes;


  /***************************************************************/
  /* print execution info */

  if (rank == 0) {
  sf_warning("==================================================");
  sf_warning("==   Lowrank two-step elastic wave propagator   ==");
  sf_warning("==================================================");
#ifdef _OPENMP
#pragma omp parallel
  {
    nth = omp_get_num_threads(); /* omp_set_num_threads(nth); */
  }
  sf_warning(">>>>>>>>>>OpenMP: Using %d threads<<<<<<<<<<", nth);
#endif
  /* wave modeling space */
  sf_warning("===== Shooting Info =====");
  if (shtnum0 != shtnum)
    sf_warning("Total shot number is not divisible by total number of nodes! shtnum padded to %d.", shtnum);
  sf_warning("shtbgn=%d shtend=%d shtint=%d shtnum=%d shtnum0=%d",shtbgn,shtend,shtint,shtnum,shtnum0);
  sf_warning("===== Modeling dimension Info =====");
  sf_warning("nx=%d nz=%d nx1=%d nz1=%d",nx,nz,nx1,nz1);
  sf_warning("x0=%f z0=%f dx=%f dz=%f",x0,z0,dx,dz);
  sf_warning("===== Receiver/Data Info =====");
  sf_warning("nt=%d dt=%f",nt,dt);
  if (mute)
    sf_warning("Muting applied: vref=%f wd=%d",vref,wd);
  else
    sf_warning("No Muting!");
  sf_warning("snpint=%d wfnt=%d wfdt=%f",snpint,wfnt,wfdt);
  sf_warning("===== Acquisition Geometry Info =====");
  sf_warning("isz=%d gpz=%d",isz,gpz);
  sf_warning("===== Source Wavelet Info =====");
  sf_warning("src=%d t0=%f f0=%f A=%f",src,t0,f0,A);
  sf_warning("===== ABC Info =====");
  sf_warning("nbt=%d nbb=%d nbl=%d nbr=%d",nbt,nbb,nbl,nbr);
  sf_warning("ct=%f cb=%f cl=%f cr=%f",ct,cb,cl,cr);
  sf_warning("===== FFT Info =====");
  if (cmplx)
    sf_warning("Using complex FFT!");
  else
    sf_warning("Using real FFT!");
  sf_warning("pad1=%d nz2=%d nx2=%d nk=%d",pad1,nz2,nx2,nk);
  sf_warning("==================================================");
  sf_warning("==   Lowrank two-step elastic wave propagator   ==");
  sf_warning("==================================================");
  }


  /***************************************************************/
  /* propagator matrices */

  if (!sf_histint(Fpl,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left for P wave",nzx);
  if (!sf_histint(Fpl,"n2",&m2p))  sf_error("Need n2= in left for P wave");
  if (!sf_histint(Fpr,"n1",&n2) || n2 != m2p) sf_error("Need n1=%d in right for P wave",m2p);
  if (!sf_histint(Fpr,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right for P wave",nk);

  if (!sf_histint(Fsl,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left for P wave",nzx);
  if (!sf_histint(Fsl,"n2",&m2s))  sf_error("Need n2= in left for P wave");
  if (!sf_histint(Fsr,"n1",&n2) || n2 != m2s) sf_error("Need n1=%d in right for P wave",m2s);
  if (!sf_histint(Fsr,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right for P wave",nk);

  plt = sf_floatalloc2(nzx,m2p);
  prt = sf_floatalloc2(m2p,nk);
  slt = sf_floatalloc2(nzx,m2s);
  srt = sf_floatalloc2(m2s,nk);  

  sf_floatread(plt[0],nzx*m2p,Fpl);
  sf_floatread(prt[0],m2p*nk,Fpr);
  sf_floatread(slt[0],nzx*m2s,Fsl);
  sf_floatread(srt[0],m2s*nk,Fsr);


  /***************************************************************/
  /* setup Output files */

  /* Data (shot record) */
  rcdz = sf_floatalloc2(nt,nx1);
  rcdx = sf_floatalloc2(nt,nx1);

  if (rank==0) {
    Fo1 = sf_output("output"); /* Elastic-wave x-component */
    Fo2 = sf_output("Recordx"); /* Elastic-wave z-component */
    puthead3x(Fo1, nt, nx1, shtnum0, dt, dx, shtint*dx, 0.0, x0, x0+shtbgn*dx);
    puthead3x(Fo2, nt, nx1, shtnum0, dt, dx, shtint*dx, 0.0, x0, x0+shtbgn*dx);

    recordz = sf_floatalloc3(nt, nx1, shtnum);
    recordx = sf_floatalloc3(nt, nx1, shtnum);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(is,ix,it)
#endif
    for (is=0; is<shtnum; is++) {
      for (ix=0; ix<nx1; ix++) {
	for (it=0; it<nt; it++) {
	  recordz[is][ix][it] = 0.0f;
	  recordx[is][ix][it] = 0.0f;
	}
      }
    }
  } else { Fo1 = NULL; Fo2 = NULL; recordz = NULL; recordx = NULL; }

  /* Wavefield Snapshot */
  wavepz = sf_floatalloc2(nzx2,wfnt);
  wavepx = sf_floatalloc2(nzx2,wfnt);
  wavesz = sf_floatalloc2(nzx2,wfnt);
  wavesx = sf_floatalloc2(nzx2,wfnt);

  if (rank==0) {
    if (NULL != sf_getstring("Pwavez")) {
      Fpz = sf_output("Pwavez");
      puthead3(Fpz, nz, nx, wfnt, dz, dx, wfdt, z0, x0, 0.0);
    } else { Fpz = NULL;}

    if (NULL != sf_getstring("Pwavex")) {
      Fpx = sf_output("Pwavex");
      puthead3(Fpx, nz, nx, wfnt, dz, dx, wfdt, z0, x0, 0.0);
    } else { Fpx = NULL; }

    if (NULL != sf_getstring("Swavez")) {
      Fsz = sf_output("Swavez");
      puthead3(Fsz, nz, nx, wfnt, dz, dx, wfdt, z0, x0, 0.0);
    } else { Fsz = NULL; }

    if (NULL != sf_getstring("Swavex")) {
      Fsx = sf_output("Swavex");
      puthead3(Fsx, nz, nx, wfnt, dz, dx, wfdt, z0, x0, 0.0);
    } else { Fsx = NULL; }
  } else { Fpz = NULL; Fpx = NULL; Fsz = NULL; Fsx = NULL; }

  if (mig) {
    /* Images */
    imgpp     = sf_floatalloc(nzx1);
    imgps     = sf_floatalloc(nzx1);
    imgsp     = sf_floatalloc(nzx1);
    imgss     = sf_floatalloc(nzx1);
    /*imgpp_sum = sf_floatalloc(nzx1);
      imgps_sum = sf_floatalloc(nzx1);
      imgsp_sum = sf_floatalloc(nzx1);
      imgss_sum = sf_floatalloc(nzx1);*/

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
    for (ix=0; ix<nzx1; ix++) {
      imgpp[ix]     = 0.0f;
      imgps[ix]     = 0.0f;
      imgsp[ix]     = 0.0f;
      imgss[ix]     = 0.0f;
      /*imgpp_sum[ix] = 0.0f;
	imgps_sum[ix] = 0.0f;
	imgsp_sum[ix] = 0.0f;
	imgss_sum[ix] = 0.0f;*/
    }

    if (rank==0) {
      if (NULL != sf_getstring("ImagePP")) {
	Fimgpp = sf_output("ImagePP");
	puthead2(Fimgpp, nz1, nx1, dz, dx, z0, x0);
      } else { Fimgpp = NULL; }

      if (NULL != sf_getstring("ImagePS")) {
	Fimgps = sf_output("ImagePS");
	puthead2(Fimgps, nz1, nx1, dz, dx, z0, x0);
      } else { Fimgps = NULL; }

      if (NULL != sf_getstring("ImageSP")) {
	Fimgsp = sf_output("ImageSP");
	puthead2(Fimgsp, nz1, nx1, dz, dx, z0, x0);
      } else { Fimgsp = NULL; }

      if (NULL != sf_getstring("ImageSS")) {
	Fimgss = sf_output("ImageSS");
	puthead2(Fimgss, nz1, nx1, dz, dx, z0, x0);
      } else { Fimgss = NULL; }
    } else { Fimgpp = NULL; Fimgps = NULL; Fimgsp = NULL; Fimgss = NULL; }
  } else {
    imgpp     = NULL;
    imgps     = NULL;
    imgsp     = NULL;
    imgss     = NULL;
    /*imgpp_sum = NULL;
      imgps_sum = NULL;
      imgsp_sum = NULL;
      imgss_sum = NULL;*/
    Fimgpp    = NULL;
    Fimgps    = NULL;
    Fimgsp    = NULL;
    Fimgss    = NULL;
  }

  /***************************************************************/
  /* set up data structure for passing info */

  geop = (geopar) sf_alloc(1, sizeof(*geop));

  geop->nx  = nx;
  geop->nz  = nz;
  geop->nx1 = nx1;
  geop->nz1 = nz1;
  geop->dx  = dx;
  geop->dz  = dz;
  geop->gpz = gpz;
/*geop->gpx = gpx;
  geop->gpl = gpl;*/
  geop->nt  = nt;
  geop->src = src;
/*geop->isx = isx;*/
  geop->isz = isz;
  geop->dt  = dt;
  geop->f0  = f0;
  geop->t0  = t0;
  geop->A   = A;
  geop->abc = abc;
  geop->nbt = nbt;
  geop->nbb = nbb;
  geop->nbl = nbl;
  geop->nbr = nbr;
  geop->ct  = ct;
  geop->cb  = cb;
  geop->cl  = cl;
  geop->cr  = cr;
  geop->m2p = m2p;
  geop->m2s = m2s;
  geop->plt = plt;
  geop->prt = prt;
  geop->slt = slt;
  geop->srt = srt;
  geop->snpint = snpint;
  geop->wfnt   = wfnt;
  geop->verb   = verb;
  geop->cmplx  = cmplx;
  geop->pad1   = pad1;
  geop->nk = nk;

  /* start to work !!! */
  tstart = clock();

  for (is=0; is*nodes<shtnum; is++){

    shtcur = is*nodes+rank; // current shot index

    if (rank==0)
      sf_warning(">>>>> Forward Propagation [ is=%d ] <<<<<",is);

    if (shtcur<shtnum0) {
      isx = shtbgn + shtint*(shtcur);
      geop->isx = isx;

      if (verb) {
	sf_warning("============================");
	sf_warning("processing shot #%d", shtcur);
	sf_warning("============================");
      }
  
      lrewefor2(rcdz, rcdx, wavepz, wavepx, wavesz, wavesx, geop);

      if (mute) {
	mutingf(nt, nx1, dt, dx, dz, isx-nbl, isz-nbt, gpz-nbt, vref, wd, rcdz);
	mutingf(nt, nx1, dt, dx, dz, isx-nbl, isz-nbt, gpz-nbt, vref, wd, rcdx);
      }
    }

    /* sum the shot gathers */
    if (rank==0) recvbuf = recordz[is*nodes][0];
    else recvbuf = NULL;
    sendbuf = rcdz[0];
    MPI_Gather(sendbuf, nx1*nt, MPI_FLOAT, recvbuf, nx1*nt, MPI_FLOAT, 0, MPI_COMM_WORLD); // recordz[is][ix][it] = rcdz[ix][it];

    if (rank==0) recvbuf = recordx[is*nodes][0];
    else recvbuf = NULL;
    sendbuf = rcdx[0];
    MPI_Gather(sendbuf, nx1*nt, MPI_FLOAT, recvbuf, nx1*nt, MPI_FLOAT, 0, MPI_COMM_WORLD); // recordx[is][ix][it] = rcdx[ix][it];

    if (shtcur==0) {
      if (NULL!=Fpz)
	for (it = 0; it < wfnt; it++)
	  for (ix = 0; ix < nx; ix++)
	    sf_floatwrite(wavepz[it]+ix*nz2,nz,Fpz);
      if (NULL!=Fpx)
	for (it = 0; it < wfnt; it++)
	  for (ix = 0; ix < nx; ix++)
	    sf_floatwrite(wavepx[it]+ix*nz2,nz,Fpx);
      if (NULL!=Fsz)
	for (it = 0; it < wfnt; it++)
	  for (ix = 0; ix < nx; ix++)
	    sf_floatwrite(wavesz[it]+ix*nz2,nz,Fsz);
      if (NULL!=Fsx)
	for (it = 0; it < wfnt; it++)
	  for (ix = 0; ix < nx; ix++)
	    sf_floatwrite(wavesx[it]+ix*nz2,nz,Fsx);
    }

    if (mig) {
      if (rank==0)
	sf_warning(">>>>> Backward Propagation [ is=%d ] <<<<<",is);

      if (shtcur<shtnum0) {
	lrewebac2(imgpp, imgps, imgsp, imgss, rcdz, rcdx, wavepz, wavepx, wavesz, wavesx, geop);
	/*
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	for (ix=0; ix<nzx1; ix++) {
	  imgpp_sum[ix] += imgpp[ix];
	  imgps_sum[ix] += imgps[ix];
	  imgsp_sum[ix] += imgsp[ix];
	  imgss_sum[ix] += imgss[ix];
	}
	*/
      }
    } /* if (mig)*/
  } /*shot iteration*/

  MPI_Barrier(MPI_COMM_WORLD);

  /******* output data *******/
  if (rank==0) {
    sf_warning("Writing Data...");
    sf_floatwrite(recordz[0][0],shtnum0*nx1*nt,Fo1);
    sf_floatwrite(recordx[0][0],shtnum0*nx1*nt,Fo2);
  }

  /******* output image *******/
  if (mig) {
    if (rank==0) {
      sf_warning("Stacking images...");
    }

    if (rank==0) {
#if MPI_VERSION >= 2
      sendbuf = (float *) MPI_IN_PLACE;
#else /* will fail */
      sendbuf = NULL;
#endif 
      recvbuf = imgpp;
    } else {
      sendbuf = imgpp;
      recvbuf = NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, nzx1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

    if (rank==0) {
#if MPI_VERSION >= 2
      sendbuf = (float *) MPI_IN_PLACE;
#else /* will fail */
      sendbuf = NULL;
#endif 
      recvbuf = imgps;
    } else {
      sendbuf = imgps;
      recvbuf = NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, nzx1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

    if (rank==0) {
#if MPI_VERSION >= 2
      sendbuf = (float *) MPI_IN_PLACE;
#else /* will fail */
      sendbuf = NULL;
#endif 
      recvbuf = imgsp;
    } else {
      sendbuf = imgsp;
      recvbuf = NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, nzx1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

    if (rank==0) {
#if MPI_VERSION >= 2
      sendbuf = (float *) MPI_IN_PLACE;
#else /* will fail */
      sendbuf = NULL;
#endif 
      recvbuf = imgss;
    } else {
      sendbuf = imgss;
      recvbuf = NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, nzx1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

    if (rank==0) {
      if (NULL!=Fimgpp) {
	sf_warning("Writing PP Image...");
	sf_floatwrite(imgpp,nzx1,Fimgpp);
      }
      if (NULL!=Fimgps) {
	sf_warning("Writing PS Image...");
	sf_floatwrite(imgps,nzx1,Fimgps);
      }
      if (NULL!=Fimgsp) {
	sf_warning("Writing SP Image...");
	sf_floatwrite(imgsp,nzx1,Fimgsp);
      }
      if (NULL!=Fimgss) {
	sf_warning("Writing SS Image...");
	sf_floatwrite(imgss,nzx1,Fimgss);
      }
    }

  } /* if (mig) */

  /*free memory*/
  free(*plt); free(plt);
  free(*prt); free(prt);
  free(*slt); free(slt);
  free(*srt); free(srt);
  free(geop);
  free(*wavepz); free(wavepz);
  free(*wavepx); free(wavepx);
  free(*wavesz); free(wavesz);
  free(*wavesx); free(wavesx);
  free(*rcdz); free(rcdz);
  free(*rcdx); free(rcdx);
  if (mig) { 
    free(imgpp); free(imgps); free(imgsp); free(imgss);
    /*free(imgpp_sum); free(imgps_sum); free(imgsp_sum); free(imgss_sum);*/
  }
  if (rank==0) {
    free(**recordz); free(*recordz); free(recordz);
    free(**recordx); free(*recordx); free(recordx);
  }

  tend = clock();
  duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
  if (verb)
    if (rank==0)
      sf_warning(">> The CPU time of single shot migration is: %f seconds << ", duration);

  MPI_Finalize();
  exit(0);
}

/*****************************************************************/
<<<<<<< .mine
=======
int fft2_init(bool cmplx1        /* if complex transform */,
	      int pad1           /* padding on the first axis */,
	      int nx,   int ny   /* input data size */, 
	      int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    /*sf_warning("Using threaded FFTW3!\n");*/
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
#else
   int i2;
#endif
	
    _cmplx = cmplx1;
	
    if (_cmplx) {
	_nk = _n1 = kiss_fft_next_fast_size(nx*pad1);
		
#ifndef SF_HAS_FFTW
	cfg1  = kiss_fft_alloc(_n1,0,NULL,NULL);
	icfg1 = kiss_fft_alloc(_n1,1,NULL,NULL);
#endif
    } else {
	_nk = kiss_fft_next_fast_size(pad1*(nx+1)/2)+1;
	_n1 = 2*(_nk-1);
		
#ifndef SF_HAS_FFTW
	cfg  = kiss_fftr_alloc(_n1,0,NULL,NULL);
	icfg = kiss_fftr_alloc(_n1,1,NULL,NULL);
#endif
    }
		
    _n2 = kiss_fft_next_fast_size(ny);

    if (_cmplx) {
	cc = sf_complexalloc2(_n1,_n2);
    } else {
	ff = sf_floatalloc2(_n1,_n2);
    }
    dd = sf_complexalloc(_nk*_n2);
	
#ifndef SF_HAS_FFTW
    cfg2  = kiss_fft_alloc(_n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(_n2,1,NULL,NULL);
 	
    tmp =    (kiss_fft_cpx **) sf_alloc(_n2,sizeof(*tmp));
    tmp[0] = (kiss_fft_cpx *)  sf_alloc(_nk*_n2,sizeof(kiss_fft_cpx));
    for (i2=0; i2 < _n2; i2++) {
	tmp[i2] = tmp[0]+i2*_nk;
    }
	
    trace2 = sf_complexalloc(_n2);
    ctrace2 = (kiss_fft_cpx *) trace2;
#endif

    *nx2 = _n1;
    *ny2 = _n2;
	
    _weight =  1.0/(_n1*_n2);
	
    return (_nk*_n2);
}

void fft2(float *inp      /* [_n1*_n2] */, 
	  sf_complex *out /* [_nk*_n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==cfg) {
	cfg = _cmplx? 
	    fftwf_plan_dft_2d(_n2,_n1,
			      (fftwf_complex *) cc[0],
			      (fftwf_complex *) dd,
			      FFTW_FORWARD, FFTW_MEASURE):
	    fftwf_plan_dft_r2c_2d(_n2,_n1,
				  ff[0], (fftwf_complex *) dd,
				  FFTW_MEASURE);
	if (NULL == cfg) sf_error("FFTW failure.");
    }
#endif

    /* FFT centering */
    for (i2=0; i2<_n2; i2++) {
	for (i1=0; i1<_n1; i1++) {
	    if (_cmplx) {
		cc[i2][i1] = sf_cmplx(((i2%2==0)==(i1%2==0))? inp[i2*_n1+i1]:-inp[i2*_n1+i1],0.);
	    } else {
		ff[i2][i1] = (i2%2)? -inp[i2*_n1+i1]:inp[i2*_n1+i1];
	    }
	}
    }
    
#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);

    for (i1=0; i1 < _nk*_n2; i1++)
      out[i1] = dd[i1];
#else	
    for (i2=0; i2 < _n2; i2++) {
	if (_cmplx) {
	    kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i2],tmp[i2],1);
	} else {
	    kiss_fftr (cfg,ff[i2],tmp[i2]);
	}
    }
	
    for (i1=0; i1 < _nk; i1++) {
	kiss_fft_stride(cfg2,tmp[0]+i1,ctrace2,_nk);
	for (i2=0; i2<_n2; i2++) {
	    out[i2*_nk+i1] = trace2[i2];
	}
    }
#endif
}

void ifft2_allocate(sf_complex *inp /* [_nk*_n2] */)
/*< allocate inverse transform >*/
{
  /* for backward compatibility */
}

void ifft2(float *out      /* [_n1*_n2] */, 
	   sf_complex *inp /* [_nk*_n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==icfg) {
      icfg = _cmplx? 
	fftwf_plan_dft_2d(_n2,_n1,
			  (fftwf_complex *) dd, 
			  (fftwf_complex *) cc[0],
			  FFTW_BACKWARD, FFTW_MEASURE):
	fftwf_plan_dft_c2r_2d(_n2,_n1,
			      (fftwf_complex *) dd, ff[0],
			      FFTW_MEASURE);
      if (NULL == icfg) sf_error("FFTW failure.");
    }
#endif

#ifdef SF_HAS_FFTW
    for (i1=0; i1 < _nk*_n2; i1++)
      dd[i1] = inp[i1];

    fftwf_execute(icfg);
#else
    for (i1=0; i1 < _nk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace2,_nk);
		
	for (i2=0; i2<_n2; i2++) {
	    tmp[i2][i1] = ctrace2[i2];
	}
    }
    for (i2=0; i2 < _n2; i2++) {
	if (_cmplx) {
	    kiss_fft_stride(icfg1,tmp[i2],(kiss_fft_cpx *) cc[i2],1);
	} else {
	    kiss_fftri(icfg,tmp[i2],ff[i2]);
	}
    }
#endif
    
    /* FFT centering and normalization */
    for (i2=0; i2<_n2; i2++) {
	for (i1=0; i1<_n1; i1++) {
	    if (_cmplx) {
		out[i2*_n1+i1] = (((i2%2==0)==(i1%2==0))? _weight:-_weight) * crealf(cc[i2][i1]);
	    } else {
		out[i2*_n1+i1] = (i2%2? -_weight: _weight)*ff[i2][i1];
	    }
	}
    }
}

void fft2_finalize()
/*< clean up fftw >*/
{
/* make sure everything is back to its pristine state */
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_cleanup_threads();
#endif
    fftwf_destroy_plan(cfg);
    fftwf_destroy_plan(icfg);
    fftwf_cleanup();
    cfg=NULL;
    icfg=NULL;
#else
    if (NULL != cfg) { free(cfg); cfg=NULL; }
    if (NULL != icfg) { free(icfg); icfg=NULL; }
    if (NULL != cfg1) { free(cfg1); cfg1=NULL; }
    if (NULL != icfg1) { free(icfg1); icfg1=NULL; }
    if (NULL != cfg2) { free(cfg2); cfg2=NULL; }
    if (NULL != icfg2) { free(icfg2); icfg2=NULL; }
    if (NULL != tmp) { free(*tmp); free(tmp); tmp=NULL; }
    if (NULL != ctrace2) { free(ctrace2); ctrace2=NULL; }
    if (NULL != trace2) { free(trace2); trace2=NULL; }
#endif
    if (NULL != ff) { free(*ff); free(ff); ff=NULL; }
    if (NULL != cc) { free(*cc); free(cc); cc=NULL; }
    if (NULL != dd) { free(dd); dd=NULL; }
}

/*****************************************************************/
void abc_init(int n1,  int n2    /*model size*/,
	      int n12, int n22   /*padded model size*/,
	      int nb1, int nb2   /*top, bottom*/,
	      int nb3, int nb4   /*left, right*/,
	      float c1, float c2 /*top, bottom*/,
	      float c3, float c4 /*left, right*/)
/*< initialization >*/
{
    int c;
    _nz = n1;
    _nx = n2;
    _nz2= n12;
    _nx2= n22;
    _nbt = nb1;
    _nbb = nb2;
    _nbl = nb3;
    _nbr = nb4;
    _ct = c1;
    _cb = c2;
    _cl = c3;
    _cr = c4;
    if(_nbt) _wt =  sf_floatalloc(_nbt);
    if(_nbb) _wb =  sf_floatalloc(_nbb);
    if(_nbl) _wl =  sf_floatalloc(_nbl);
    if(_nbr) _wr =  sf_floatalloc(_nbr);
    c=0;
    abc_cal(c,_nbt,_ct,_wt);
    abc_cal(c,_nbb,_cb,_wb);
    abc_cal(c,_nbl,_cl,_wl);
    abc_cal(c,_nbr,_cr,_wr);
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(_nbt) free(_wt);
    if(_nbb) free(_wb);
    if(_nbl) free(_wl);
    if(_nbr) free(_wr);
}

void abc_apply(float *a /*2-D matrix*/) 
/*< boundary decay>*/
{
    int i;
    int iz, ix;

    /* top */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix,i)
{
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbt; iz++) {  
        for (ix=_nbl; ix < _nx-_nbr; ix++) {
	  i = _nz2*ix + iz;
	  a[i] *= _wt[iz];
        }
    }
    /* bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbb; iz++) {  
        for (ix=_nbl; ix < _nx-_nbr; ix++) {
	  i = _nz2*ix + _nz-1-iz;
	  a[i] *= _wb[iz];
        }
    }
    /* left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=_nbt; iz < _nz-_nbb; iz++) {  
        for (ix=0; ix < _nbl; ix++) {
	  i = _nz2*ix + iz;
	  a[i] *= _wl[ix];
        }
    }
    /* right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=_nbt; iz < _nz-_nbb; iz++) {  
        for (ix=0; ix < _nbr; ix++) {
	  i = _nz2*(_nx-1-ix) + iz;
          a[i] *= _wr[ix];
        }
    }
    /* top left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbt; iz++) {  
        for (ix=0; ix < _nbl; ix++) {
	  i = _nz2*ix + iz;
	  a[i] *= iz>ix? _wl[ix]:_wt[iz];
        }
    }
    /* top right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbt; iz++) {  
        for (ix=0; ix < _nbr; ix++) {
	  i = _nz2*(_nx-1-ix) + iz;
	  a[i] *= iz>ix? _wr[ix]:_wt[iz];
        }
    }
    /* bottom left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbb; iz++) {  
        for (ix=0; ix < _nbl; ix++) {
	  i = _nz2*ix + _nz-1-iz;
          a[i] *= iz>ix? _wl[ix]:_wb[iz];
        }
    }
    /* bottom right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < _nbb; iz++) {  
        for (ix=0; ix < _nbr; ix++) {
	  i = _nz2*(_nx-1-ix) + _nz-1-iz;
          a[i] *= iz>ix? _wr[ix]:_wb[iz];
        }
    }

#ifdef _OPENMP
}
#endif
}

void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    /*const float pi=SF_PI;*/
    if(!nb) return;
    switch(abc) {
    default:
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
        for(ib=0; ib<nb; ib++){
	    w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	}
    }
}

/*****************************************************************/
void puthead3(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (z,x,t) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","t");
        sf_putstring(Fo,"unit1","m");
        sf_putstring(Fo,"unit2","m");
        sf_putstring(Fo,"unit3","second");
}

void puthead3x(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (t,x,s) shot gather cubes  >*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","t");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","shot");
        sf_putstring(Fo,"unit1","second");
        sf_putstring(Fo,"unit2","m");
        sf_putstring(Fo,"unit3","m");
}

void puthead2(sf_file Fo, int n1, int n2, float d1, float d2, float o1, float o2)
/*<  put head for (x,z) domain float-type 2D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","m");
        sf_putstring(Fo,"unit2","m");
}

void puthead2x(sf_file Fo, int n1, int n2, float d1, float d2, float o1, float o2)
/*<  put head for (x,z) domain float-type 2D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","t");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","second");
        sf_putstring(Fo,"unit2","m");
}

/*****************************************************************/
float Ricker(float t, float f0, float t0, float A) 
/*< ricker wavelet:
 * f0: peak frequency
 * t0: time lag
 * A: amplitude
 * ************************>*/
{
        float x=pow(SF_PI*f0*(t-t0),2);
        return -A*exp(-x)*(1-2*x);
}

/*****************************************************************/
>>>>>>> .r13904
/* forward modeling operator: generating shot gather */
int lrewefor2(float **rcdz, float** rcdx, float **upz, float **upx, float **usz, float **usx, geopar geop)
/*< lowrank 2d elastic forward modeling >*/
{
  int nx,nz,nz1,nx1;     /* model dimension parameters */
  float dx,dz;           /* spatial sampling */
  int gpz;

  int nt, src, isx, isz; /* model/data dimension and source location */
  float dt, f0, t0, A;   /* time and source related */

  bool abc;              /* abc flag */
  int nbt,nbb,nbl,nbr;   /* abc width */
  float ct,cb,cl,cr;     /* decay parameter */

  int m2p,m2s;           /* lowrank matrices dimension */
  float **plt, **prt, **slt, **srt; /* lowrank matrices */

  int snpint, wfnt; /* snapshot interval and total number */

  bool verb;  /* verbose flag */
  bool cmplx; /* complex fft flag */
  int pad1;   /* fft padding on the first axis */

  /* computation related variables/arrays*/
  float t;
  float c, old;
  int wfit;
  int iz, ix, ik, im, it, i, j;
  int nzx,nzx2,nz2,nx2,nk,nkz,nkx;
  float kz0, kx0, dkz, dkx;
  float *kx,*kz,kmod;
  float **wavep, **waves, *uz2, *uz1, *ux2, *ux1;
  sf_complex *cwavepz, *cwavepx, *cwavesz, *cwavesx, *cwavex, *cwavez, *cwavem;

  nx  = geop->nx;
  nz  = geop->nz;
  nx1 = geop->nx1;
  nz1 = geop->nz1;
  dx  = geop->dx;
  dz  = geop->dz;
  gpz = geop->gpz;
  //gpx = geop->gpx;
  //gpl = geop->gpl;
  nt  = geop->nt;
  src = geop->src;
  isx = geop->isx;
  isz = geop->isz;
  dt  = geop->dt;
  f0  = geop->f0;
  t0  = geop->t0;
  A   = geop->A;
  abc = geop->abc;
  nbt = geop->nbt;
  nbb = geop->nbb;
  nbl = geop->nbl;
  nbr = geop->nbr;
  ct  = geop->ct;
  cb  = geop->cb;
  cl  = geop->cl;
  cr  = geop->cr;
  m2p = geop->m2p;
  m2s = geop->m2s;
  plt = geop->plt;
  prt = geop->prt;
  slt = geop->slt;
  srt = geop->srt;
  snpint = geop->snpint;
  wfnt   = geop->wfnt;
  verb   = geop->verb;
  cmplx  = geop->cmplx;
  pad1   = geop->pad1;

  nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
  if (nk!=geop->nk) sf_error("Incompatible nk!");
  nzx = nz*nx;
  nzx2 = nz2*nx2;
  kz0 = cmplx ? -SF_PI/dz : 0;
  kx0 = -SF_PI/dx;
  dkz = 2*SF_PI/(dz*nz2);
  dkx = 2*SF_PI/(dx*nx2);
  nkz = cmplx ? nz2 : 1+nz2/2;
  nkx = nx2;

  kx=sf_floatalloc(nk);
  kz=sf_floatalloc(nk);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik,kmod)
#endif
  for (ix=0; ix < nkx; ix++) {
    for (iz=0; iz < nkz; iz++) {
      ik = iz + ix*nkz;
      kz[ik] = kz0+iz*dkz;
      kx[ik] = kx0+ix*dkx;
      kmod = sqrtf(kz[ik]*kz[ik]+kx[ik]*kx[ik]);
      kz[ik] /= kmod;
      kx[ik] /= kmod;
    }
  }

  /**************** wavefield propagation ****************/

  uz1 = sf_floatalloc(nzx2);
  ux1 = sf_floatalloc(nzx2);
  uz2 = sf_floatalloc(nzx2);
  ux2 = sf_floatalloc(nzx2);

  cwavepz= sf_complexalloc(nk);
  cwavepx= sf_complexalloc(nk);
  cwavesz= sf_complexalloc(nk);
  cwavesx= sf_complexalloc(nk);
  cwavez = sf_complexalloc(nk);
  cwavex = sf_complexalloc(nk);
  cwavem = sf_complexalloc(nk);
  wavep  = sf_floatalloc2(nzx2,m2p);
  waves  = sf_floatalloc2(nzx2,m2s);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
  for (iz=0; iz < nzx2; iz++) {
    uz1[iz]=0.;
    ux1[iz]=0.;
    uz2[iz]=0.;
    ux2[iz]=0.;
  }

  /* initialize abc */
  if (abc)
    abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

  //  ifft2_allocate(cwavem);

  wfit = 0;
  for(it=0;it<nt;it++) {
    t=it*dt;

    /* forward-propagating using two-step k-space solution in VTI media */

    /* Source Injection */
    if (src==1) {
      // 2D exploding force source
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
      for(ix=-1;ix<=1;ix++) {
	for(iz=-1;iz<=1;iz++) {
	  if(fabs(ix)+fabs(iz)==2) {
	    j = isz+iz+nz2*(isx+ix);
	    uz2[j]+=iz*Ricker(t, f0, t0, A);
	    ux2[j]+=ix*Ricker(t, f0, t0, A);
	  }
	}
      }
    } else {
      // 2D equil-energy force source
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
      for(ix=-1;ix<=1;ix++) {
        for(iz=-1;iz<=1;iz++) {
	  if(fabs(ix)+fabs(iz)==2) {
	    j = isz+iz+nz2*(isx+ix);
	    if(ix==-1&&iz==1)  
	      uz2[j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==-1&&iz==-1) 
	      ux2[j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==1&&iz==1)  
	      ux2[j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==1&&iz==-1) 
	      uz2[j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	  }
        }
      }
    }

    /* Going to wavenumber domain */
    fft2(uz2,cwavez);
    fft2(ux2,cwavex);

    /* First decompose the wavefield into P and S wave components */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik)
#endif
    for (ix=0; ix < nkx; ix++) {
      for (iz=0; iz < nkz; iz++) {
	ik = iz + ix*nkz;
	cwavepz[ik] = cwavez[ik]*kz[ik]*kz[ik] + cwavex[ik]*kx[ik]*kz[ik]; /* upz(k) */
	cwavesz[ik] = cwavez[ik] - cwavepz[ik]; /* usz(k) */
	cwavepx[ik] = cwavez[ik]*kz[ik]*kx[ik] + cwavex[ik]*kx[ik]*kx[ik]; /* upx(k) */
	cwavesx[ik] = cwavex[ik] - cwavepx[ik]; /* usx(k) */
      }
    }

    /* writing to wavefield "cube" */
    if (it%snpint == 0) {
      if (NULL != upz)
        ifft2(upz[wfit],cwavepz);
      if (NULL != upx)
        ifft2(upx[wfit],cwavepx);
      if (NULL != usz)
        ifft2(usz[wfit],cwavesz);
      if (NULL != usx)
        ifft2(usx[wfit],cwavesx);
      wfit ++;
    }

    /* Computing the Z component */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepz[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepz[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesz[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesz[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = uz2[j];
	c += c - uz1[j];
	uz1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}

	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	uz2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(uz1);
      abc_apply(uz2);
    }

    /* Computing the X component */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepx[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepx[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesx[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesx[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = ux2[j];
	c += c - ux1[j];
	ux1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}
	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	ux2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(ux1);
      abc_apply(ux2);
    }

    /* record data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,j)
#endif
    for (ix = 0; ix < nx1; ix++) {
      j = (ix+nbl)*nz2 + gpz;
      rcdz[ix][it] = uz2[j];
      rcdx[ix][it] = ux2[j];
    }

    if(it%100==0)
      if (verb)
	sf_warning("Elastic: it= %d",it);
  }/* it loop */

  if (wfit!=wfnt) sf_warning("wfnt discrepancy! wfit = %d",wfit);

  /* clean up */
  
  if (abc)
    abc_close();

  fft2_finalize();

  free(kx); free(kz);
  free(*wavep); free(wavep); free(*waves); free(waves);
  free(uz2); free(uz1); free(ux2); free(ux1);
  free(cwavepz); free(cwavepx); free(cwavesz); free(cwavesx);
  free(cwavex); free(cwavez);
  free(cwavem);

  return 0;
}

/* backward propagation and applying imaging condition: generate image */
int lrewebac2(float *imgpp, float *imgps, float *imgsp, float *imgss, float **rcdz, float** rcdx, float **upz, float **upx, float **usz, float **usx, geopar geop)
/*< lowrank 2d elastic backward modeling >*/
{
  int nx,nz,nz1,nx1;     /* model dimension parameters */
  float dx,dz;           /* spatial sampling */
  int gpz;

  int nt, src, isx, isz; /* model/data dimension and source location */
  float dt, f0, t0, A;   /* time and source related */

  bool abc;              /* abc flag */
  int nbt,nbb,nbl,nbr;   /* abc width */
  float ct,cb,cl,cr;     /* decay parameter */

  int m2p,m2s;           /* lowrank matrices dimension */
  float **plt, **prt, **slt, **srt; /* lowrank matrices */

  int snpint, wfnt; /* snapshot interval and total number */

  bool verb;  /* verbose flag */
  bool cmplx; /* complex fft flag */
  int pad1;   /* fft padding on the first axis */

  /* computation related variables/arrays*/
  float t;
  float c, old;
  int wfit;
  int iz, ix, ik, im, it, i, j;
  int nzx,nzx2,nz2,nx2,nk,nkz,nkx;
  float kz0, kx0, dkz, dkx;
  float *kx,*kz,kmod;
  float **wavep, **waves, *uz2, *uz1, *ux2, *ux1;
  sf_complex *cwavepz, *cwavepx, *cwavesz, *cwavesx, *cwavex, *cwavez, *cwavem;
  float *utmpz, *utmpx; /* temporary wavefield for imaging condition*/

  nx  = geop->nx;
  nz  = geop->nz;
  nx1 = geop->nx1;
  nz1 = geop->nz1;
  dx  = geop->dx;
  dz  = geop->dz;
  gpz = geop->gpz;
  //gpx = geop->gpx;
  //gpl = geop->gpl;
  nt  = geop->nt;
  src = geop->src;
  isx = geop->isx;
  isz = geop->isz;
  dt  = geop->dt;
  f0  = geop->f0;
  t0  = geop->t0;
  A   = geop->A;
  abc = geop->abc;
  nbt = geop->nbt;
  nbb = geop->nbb;
  nbl = geop->nbl;
  nbr = geop->nbr;
  ct  = geop->ct;
  cb  = geop->cb;
  cl  = geop->cl;
  cr  = geop->cr;
  m2p = geop->m2p;
  m2s = geop->m2s;
  plt = geop->plt;
  prt = geop->prt;
  slt = geop->slt;
  srt = geop->srt;
  snpint = geop->snpint;
  wfnt   = geop->wfnt;
  verb   = geop->verb;
  cmplx  = geop->cmplx;
  pad1   = geop->pad1;

  nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
  if (nk!=geop->nk) sf_error("Incompatible nk!");
  nzx = nz*nx;
  nzx2 = nz2*nx2;
  kz0 = cmplx ? -SF_PI/dz : 0;
  kx0 = -SF_PI/dx;
  dkz = 2*SF_PI/(dz*nz2);
  dkx = 2*SF_PI/(dx*nx2);
  nkz = cmplx ? nz2 : 1+nz2/2;
  nkx = nx2;

  kx=sf_floatalloc(nk);
  kz=sf_floatalloc(nk);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik,kmod)
#endif
  for (ix=0; ix < nkx; ix++) {
    for (iz=0; iz < nkz; iz++) {
      ik = iz + ix*nkz;
      kz[ik] = kz0+iz*dkz;
      kx[ik] = kx0+ix*dkx;
      kmod = sqrtf(kz[ik]*kz[ik]+kx[ik]*kx[ik]);
      kz[ik] /= kmod;
      kx[ik] /= kmod;
    }
  }

  /**************** wavefield propagation ****************/

  uz1 = sf_floatalloc(nzx2);
  ux1 = sf_floatalloc(nzx2);
  uz2 = sf_floatalloc(nzx2);
  ux2 = sf_floatalloc(nzx2);
  utmpz = sf_floatalloc(nzx2);
  utmpx = sf_floatalloc(nzx2);

  cwavepz= sf_complexalloc(nk);
  cwavepx= sf_complexalloc(nk);
  cwavesz= sf_complexalloc(nk);
  cwavesx= sf_complexalloc(nk);
  cwavez = sf_complexalloc(nk);
  cwavex = sf_complexalloc(nk);
  cwavem = sf_complexalloc(nk);
  wavep  = sf_floatalloc2(nzx2,m2p);
  waves  = sf_floatalloc2(nzx2,m2s);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
  for (iz=0; iz < nzx2; iz++) {
    uz1[iz]=0.;
    ux1[iz]=0.;
    uz2[iz]=0.;
    ux2[iz]=0.;
    utmpz[iz]=0.;
    utmpx[iz]=0.;
  }

  /* initialize abc */
  if (abc)
    abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

  //  ifft2_allocate(cwavem);

  wfit = wfnt-1;
  for(it=nt-1;it>=0;it--) {
    t=it*dt;

    /* backward-propagating using two-step k-space solution in VTI media */

    /* Data Injection */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,j)
#endif
    for (ix = 0; ix < nx1; ix++) {
      j = (ix+nbl)*nz2 + gpz;
      uz2[j] += rcdz[ix][it];
      ux2[j] += rcdx[ix][it];
    }

    /* Going to wavenumber domain */
    fft2(uz2,cwavez);
    fft2(ux2,cwavex);

    /* First decompose the wavefield into P and S wave components */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik)
#endif
    for (ix=0; ix < nkx; ix++) {
      for (iz=0; iz < nkz; iz++) {
	ik = iz + ix*nkz;
	cwavepz[ik] = cwavez[ik]*kz[ik]*kz[ik] + cwavex[ik]*kx[ik]*kz[ik]; /* upz(k) */
	cwavesz[ik] = cwavez[ik] - cwavepz[ik]; /* usz(k) */
	cwavepx[ik] = cwavez[ik]*kz[ik]*kx[ik] + cwavex[ik]*kx[ik]*kx[ik]; /* upx(k) */
	cwavesx[ik] = cwavex[ik] - cwavepx[ik]; /* usx(k) */
      }
    }

    /* applying imaging condition */
    if (it%snpint == 0) {
      if (0) {
      if (NULL != upz)
        ifft2(upz[wfit],cwavepz);
      if (NULL != upx)
        ifft2(upx[wfit],cwavepx);
      if (NULL != usz)
        ifft2(usz[wfit],cwavesz);
      if (NULL != usx)
        ifft2(usx[wfit],cwavesx);
      wfit --;
      } else {
      /* XP imaging condition: using recerver P wavefield */
      ifft2(utmpz,cwavepz);
      ifft2(utmpx,cwavepx);
      if (NULL != imgpp) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
        for (ix = 0; ix < nx1; ix++) {
          for (iz = 0; iz < nz1; iz++) {
            i = ix*nz1 + iz;
            j = (ix+nbl)*nz2 + iz+nbt;
            imgpp[i] += upz[wfit][j]*utmpz[j] + upx[wfit][j]*utmpx[j];
          }
        }
      }
      if (NULL != imgsp) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
        for (ix = 0; ix < nx1; ix++) {
          for (iz = 0; iz < nz1; iz++) {
            i = ix*nz1 + iz;
            j = (ix+nbl)*nz2 + iz+nbt;
            imgsp[i] += usz[wfit][j]*utmpz[j] + usx[wfit][j]*utmpx[j];
          }
        }
      }
      /* XS imaging condition: using recerver S wavefield */
      ifft2(utmpz,cwavesz);
      ifft2(utmpx,cwavesx);
      if (NULL != imgps) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
        for (ix = 0; ix < nx1; ix++) {
          for (iz = 0; iz < nz1; iz++) {
            i = ix*nz1 + iz;
            j = (ix+nbl)*nz2 + iz+nbt;
            imgps[i] += upz[wfit][j]*utmpz[j] + upx[wfit][j]*utmpx[j];
          }
        }
      }
      if (NULL != imgss) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
        for (ix = 0; ix < nx1; ix++) {
          for (iz = 0; iz < nz1; iz++) {
            i = ix*nz1 + iz;
            j = (ix+nbl)*nz2 + iz+nbt;
            imgss[i] += usz[wfit][j]*utmpz[j] + usx[wfit][j]*utmpx[j];
          }
        }
      }
      wfit --;
      }
    }

    /* Computing the Z component */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepz[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepz[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesz[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesz[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = uz2[j];
	c += c - uz1[j];
	uz1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}

	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	uz2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(uz1);
      abc_apply(uz2);
    }

    /* Computing the X component */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepx[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepx[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesx[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesx[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = ux2[j];
	c += c - ux1[j];
	ux1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}
	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	ux2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(ux1);
      abc_apply(ux2);
    }

    if(it%100==0)
      if (verb)
	sf_warning("Elastic: it= %d",it);
  }/* it loop */

  if (wfit!=-1) sf_warning("wfnt discrepancy! wfit = %d",wfit);

  /* clean up */
  
  if (abc)
    abc_close();

  fft2_finalize();

  free(kx); free(kz);
  free(*wavep); free(wavep); free(*waves); free(waves);
  free(uz2); free(uz1); free(ux2); free(ux1);
  free(cwavepz); free(cwavepx); free(cwavesz); free(cwavesx);
  free(cwavex); free(cwavez);
  free(cwavem);
  free(utmpz); free(utmpx);

  return 0;
}

