/* isotropic reverse-time migration */
/*
  Copyright (C) 2012 KAUST

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

/* ----------------------------------- */
/* Author : Xuxin Ma 2012              */
/* ----------------------------------- */

/* ------------ Equations ----------- */
/*                                    */
/* dp/dt = v^2 * (du/dx + dw/dz + f)  */
/* du/dt = dp/dx                      */
/* dw/dt = dp/dz                      */
/*                                    */
/* p     : stress                     */
/* u,w   : particle momentum * (-1)   */
/* v     : P-wave velocity            */
/* f  : volume source                 */
/* ---------------------------------- */

/* ------------------------------------------------------------ */
/* I/O        File      Type   Size          Description        */
/* ------------------------------------------------------------ */
/* [In/Out]   par.txt   ascii                parameters         */
/* [In]       velo.bin  float      [nx][nz]  velocity           */
/* [In]       sour.bin  float          [nt]  source function    */
/* [In/Out]   data.bin  float  [ns][nh][nt]  shot gathers       */
/* [In/Out]   uus.bin   float  [n3][Nx][Nz]  source wavefield   */
/* [Out]      uur.bin   float  [n3][Nx][Nz]  receiver wavefield */
/* [Out]      img.bin   float      [Nx][Nz]  image              */
/* [Out]      cig.bin   float  [nH][nC][Nz]  offset cig         */
/*                                                              */
/* * mode=0 : both                                              */
/*   mode=1 : source wavefield only                             */
/*   mode=2 : receiver wavefield only                           */
/*   mode=3 : source wavefield only (generate dataset)          */
/* * data.bin not needed if mode=1 or 3                         */
/* * sour.bin not needed if mode=2                              */
/* * need round (nt-1) / (n3-1)                                 */
/* * stability requirement :  dx/dt >= 5 * vmax                 */
/* * one shot per MPI process                                   */
/* * sponge boundary (typical b=40 c=2e-5)                      */
/* * cross-correlation imaging condition                        */
/* ------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <fftw3.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <rsf.h>

#include "mpiutil.h"
#include "mpisponge.h"
#include "interp.h"

#undef DEBUG

#define PREFIX_UUS "uus"
#define PREFIX_UUR "uur"
#define PREFIX_IMG "img"
#define PREFIX_CIG "cig"

#define MASTER 0
#define TAG_CS 1
#define TAG_CR 2
#define TAG_DA 4
#define TAG_NR 8

int mpi_size,mpi_rank;
int omp_nth;

int nz; float dz,z0,*z,*kz; int bzl,bzh,Nz;
int nx; float dx,x0,*x,*kx; int bxl,bxh,Nx;
int nt; float dt,   *t;
int n3; float d3;           int j3;
int nH; int   dH;    /* subsurface offset */
int nC; int   dC,C0; /* CIG position      */
int nr; float *cr;   /* [nr][2] receiver  */
float  cs[2];/* [1 ][2] shot      */

int nr,nxz,Nxz,nrt,NHCz;

Sponge *sponge;

/* char fname_par []="par.txt"; */
/* char fname_velo[]="velo.bin"; */
sf_file fvelo; float *velo; 
/* char fname_sour[]="sour.bin"; */
sf_file fsour; float *sour;
char fname_data[]="data.bin"; float *data,*datat;
char fname_uus[128];          float *uus;
char fname_uur[128];
char fname_img[128];          float *img;
char fname_cig[128];          float *cig;

int forward;
int backward;
int savedata;

/* ======== tau domain parameters ======== */
int tau; /* if true, tau domain */
char fname_vmap[]="vmap.bin"; float *vmap;
char fname_sigm[]="sigm.bin"; float *sigm;

static void init(int argc, char* argv[])
{
    int i,j,ns,nh,b[4],n[3],mode;
    float czl,czh,cxl,cxh,s0,ds,h0,dh,zs=0.0f,zr=0.0f,z_,x_,cs_w[2],*cr_w,c[4],*temp;
    FILE *Fdata;
    MPI_Status status;

    /* default values */
    nz = 1; dz = 0; z0 = 0;
    nx = 1; dx = 0; x0 = 0;
    nt = 1; dt = 0; j3 = 1;
    nh = 1; h0 = 0; dh = 0;
    ns = 1; s0 = 0; ds = 0;
    nH = 1;         dH = 1;
    nC = 1; C0 = 1; dC = 1;
    bzl = bzh = bxl = bxh = 40;
    czl = czh = cxl = cxh = 2.e-5;  
    mode = 0; tau = 0;

    if (MASTER == mpi_rank) {
	sf_init(argc,argv);

	sf_getint("nz",&nz);
	sf_getint("nz",&nx);
	sf_getint("tz",&nt);

	sf_getfloat("dz",&dz);
	sf_getfloat("dz",&dx);
	sf_getfloat("dt",&dt);

	sf_getfloat("z0",&z0);
	sf_getfloat("x0",&x0);

	sf_getint("bzl",&bzl);
	sf_getint("bzh",&bzl);
	sf_getint("bxl",&bzl);
	sf_getint("bxh",&bzl);
	
	sf_getfloat("czl",&czl);
	sf_getfloat("czh",&czl);
	sf_getfloat("cxl",&czl);
	sf_getfloat("cxh",&czl);

	sf_getint("j3",&j3);

	sf_getint("nh",&nh);
	sf_getint("ns",&ns);

	sf_getfloat("dh",&dh);
	sf_getfloat("ds",&ds);

	sf_getfloat("h0",&h0);
	sf_getfloat("s0",&s0);

	sf_getfloat("zr",&zr);
	sf_getfloat("zs",&zs);

	sf_getint("nH",&nH);
	sf_getint("nC",&nC);

	sf_getint("dH",&dH);
	sf_getint("dC",&dC);
	sf_getint("C0",&C0);

	sf_getint("mode",&mode);
	sf_getint("tau",&tau);

        if (mpi_size != ns) exit_error("Need one shot per process");
    }

    /* set options */

    if (MASTER == mpi_rank) {
	switch (mode) {
	    case 0 : forward=1; backward=1; savedata=0; break;
	    case 1 : forward=1; backward=0; savedata=0; break;
	    case 2 : forward=0; backward=1; savedata=0; break;
	    case 3 : forward=1; backward=0; savedata=1; break;
	    default: exit_error("unknown mode=%d",mode);
	}
    }

    MPI_Bcast(&forward ,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&backward,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&savedata,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&tau     ,1,MPI_INT,MASTER,MPI_COMM_WORLD);

    /* set basic dimensions */

    MPI_Bcast(&nz,1,MPI_INT  ,MASTER,MPI_COMM_WORLD); 
    MPI_Bcast(&nx,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&nt,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&dz,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&dx,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&dt,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);

    /* update boundaries for faster FFT */

    if (MASTER == mpi_rank) {
        fft_expand(nz,&bzl,&bzh);
        fft_expand(nx,&bxl,&bxh);
    }

    MPI_Bcast(&bzl,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&bzh,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&bxl,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&bxh,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&czl,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&czh,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&cxl,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&cxh,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);

    Nz = nz + bzl + bzh;
    Nx = nx + bxl + bxh;

    nxz = nx * nz;
    Nxz = Nx * Nz;

    /* check cig positions */

    if (MASTER == mpi_rank) {
        if (C0 < 0   || C0 + dC*(nC-1) < 0 ||
	    C0 >= nx || C0 + dC*(nC-1) >= nx)
	    exit_error("need cig position inside domain");
    }

    MPI_Bcast(&nH,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&dH,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&nC,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&C0,1,MPI_INT,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&dC,1,MPI_INT,MASTER,MPI_COMM_WORLD);

    NHCz = nH * nC * Nz;

    /* spatial axis */

    z = (float *)calloc(Nz,sizeof(float));
    x = (float *)calloc(Nx,sizeof(float));

    for (i=0; i < Nz; i++) z[i] = z0 + dz * (i-bzl);
    for (i=0; i < Nx; i++) x[i] = x0 + dx * (i-bxl);

    /* wavenumber axis */

    kz = (float *)calloc(Nz,sizeof(float));
    kx = (float *)calloc(Nx,sizeof(float));

    compute_k(kz,Nz);
    compute_k(kx,Nx);

    /* wavefield snapshot axis */

    if (MASTER == mpi_rank) {
        if ((nt-1) % j3) printf("** (nt-1) / j3 not round  **\n");
        n3 = 1 + (nt-1) / j3;
        d3 = dt * j3;
    }

    MPI_Bcast(&j3,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&n3,1,MPI_INT  ,MASTER,MPI_COMM_WORLD);
    MPI_Bcast(&d3,1,MPI_FLOAT,MASTER,MPI_COMM_WORLD);

    /* write new parameters 

    if (MASTER == mpi_rank) {
        fseek(Fpar,0,SEEK_END);
        fprintf(Fpar,"out : Nz=%d\n", Nz);
        fprintf(Fpar,"out : Nx=%d\n", Nx);
        fprintf(Fpar,"out : bzl=%d\n",bzl);
        fprintf(Fpar,"out : bzh=%d\n",bzh);
        fprintf(Fpar,"out : bxl=%d\n",bxl);
        fprintf(Fpar,"out : bxh=%d\n",bxh);
        fprintf(Fpar,"out : z[ 0]=%g\n", z[0]);
        fprintf(Fpar,"out : x[ 0]=%g\n", x[0]);
        fprintf(Fpar,"out : z[-1]=%g\n", z[Nz-1]);
        fprintf(Fpar,"out : x[-1]=%g\n", x[Nx-1]);
        fprintf(Fpar,"out : n3=%d\n" ,n3);
        fprintf(Fpar,"out : d3=%g\n" ,d3);
        fclose(Fpar);
    }

    */

    /* initialize sponge boundary */

    b[0] = bzl; b[1] = bzh; b[2] = bxl; b[3] = bxh;
    c[0] = czl; c[1] = czh; c[2] = cxl; c[3] = cxh;
    sponge = sponge_init(b,c);

    /* read velocity */

    velo = sf_floatalloc(nxz);

    if (MASTER == mpi_rank) {
	fvelo = sf_input("velo");
        sf_floatread(velo,nxz,fvelo);
	sf_fileclose(fvelo);
    }

    MPI_Bcast(velo,nxz,MPI_FLOAT,MASTER,MPI_COMM_WORLD);

    /* read vmap and sigm */
    if (tau) {
        vmap = (float *)calloc(nxz,sizeof(float));
        sigm = (float *)calloc(nxz,sizeof(float));

        if (MASTER == mpi_rank) {
            read_float(vmap,nxz,fname_vmap);
            read_float(sigm,nxz,fname_sigm);
        }

        MPI_Bcast(vmap,nxz,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
        MPI_Bcast(sigm,nxz,MPI_FLOAT,MASTER,MPI_COMM_WORLD);
    }

    sprintf(fname_uus,"%s%d.bin",PREFIX_UUS,mpi_rank);
    sprintf(fname_uur,"%s%d.bin",PREFIX_UUR,mpi_rank);
    sprintf(fname_img,"%s%d.bin",PREFIX_IMG,mpi_rank);
    sprintf(fname_cig,"%s%d.bin",PREFIX_CIG,mpi_rank);

#ifdef _OPENMP
#pragma omp parallel
    omp_nth = omp_get_num_threads();
    printf("process %3d : omp using %d threads\n",mpi_rank,omp_nth);
    printf("process %3d : forward=%d backward=%d savedata=%d\n",mpi_rank,forward,backward,savedata);
#endif

    /* -------------------------------- */
    /* initialize forward extrapolation */
    /* -------------------------------- */

    if (forward) {

        /* shot position */

        if (MASTER == mpi_rank) {
            cs_w[0] = cs[0] = zs;        // z
            cs_w[1] = (cs[1] = s0) + ds; // x

            for (i=1; i < mpi_size; i++, cs_w[1] += ds)
		MPI_Send(cs_w,2,MPI_FLOAT,i,TAG_CS,MPI_COMM_WORLD);
        } else {
            MPI_Recv(cs,2,MPI_FLOAT,MASTER,TAG_CS,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status,MPI_FLOAT,&j);
            if (2 != j) exit_error("failed receiving cs. %d expected %d received",2,j);
        }

        z_ = cs[0];
        x_ = cs[1];
	if (z_ < z[0] || z_ > z[Nz-1] ||
            x_ < x[0] || x_ > x[Nx-1])
	    exit_error("shot %d at x=%f z=%f outside domain",i,x_,z_);

        /* source signature */

        sour = sf_floatalloc(nt);
        if (MASTER == mpi_rank) {
	    fsour = sf_input("sour");
	    sf_floatread(sour,nt,fsour);
	    sf_fileclose(fsour);
	}
        MPI_Bcast(sour,nt,MPI_FLOAT,MASTER,MPI_COMM_WORLD);

        printf("process %3d : shot at x=%8.2f z=%8.2f\n",mpi_rank,x_,z_);
    }

    /* --------------------------------- */
    /* initialize backward extrapolation */
    /* --------------------------------- */

    /* read -> cut -> transpose */
    /* shall be re-organized. for example, MASTER cut each shot gather then send to workers */

    /* offset to receiver positions */

    if (backward || savedata) {
        if (MASTER == mpi_rank) nr = nh;
        MPI_Bcast(&nr,1,MPI_INT,MASTER,MPI_COMM_WORLD);

        cr = (float *)calloc(nr*2,sizeof(float));
        if (MASTER == mpi_rank) {
            for (i=0; i < nr; i++) {
                cr[i*2    ] = zr;
                cr[i*2 + 1] = s0 + h0 + i*dh;
            }

            cr_w = (float *)calloc(nr*2,sizeof(int));
            for (i=0; i < nr*2; i++) cr_w[i] = cr[i];

            for (i=1; i < mpi_size; i++) {
                for (j=0; j < nr; j++) cr_w[j*2+1] += ds;
                MPI_Send(cr_w,nr*2,MPI_FLOAT,i,TAG_CR,MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(cr,nr*2,MPI_FLOAT,MASTER,TAG_CR,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status,MPI_INT,&j);
            if (nr*2 != j) exit_error("failed receiving cr. %d expected %d received",nr*2,j);
        }
    }

    /* shot gathers */

    if (backward) {
        datat= (float *)calloc(nrt = nr * nt,sizeof(float)); // nr uncutted

        if (MASTER == mpi_rank) {
            if (NULL == (Fdata = fopen(fname_data,"r")))
                exit_error("failed to open %s",fname_data);

            if (nrt != fread(datat,sizeof(float),nrt,Fdata))
                exit_error("failed reading %s",fname_data);           

	    temp = (float *)calloc(nrt,sizeof(float));
            for (i=1; i < mpi_size; i++) {
                if (nrt != fread(temp,sizeof(float),nrt,Fdata))
                    exit_error("failed reading %s",fname_data);
                if (MPI_SUCCESS != MPI_Send(temp,nrt,MPI_FLOAT,i,TAG_DA,MPI_COMM_WORLD))
		    exit_error("failed sending data");
            }
            free(temp);
        } else {
            if (MPI_SUCCESS != MPI_Recv(datat,nrt,MPI_FLOAT,MASTER,TAG_DA,MPI_COMM_WORLD,&status))
		exit_error("failed receiving data");
            MPI_Get_count(&status,MPI_FLOAT,&j);
            if (nrt != j) exit_error("failed receiving data. %d expected %d received",nrt,j);
        }
    }

    if (backward || savedata) {
        n[0] = nr; n[1] = nt;
        c[0] = x0; c[1] = x0 + (nx-1)*dx;
        nr = cut_offset(backward ? datat : NULL,cr,n,c);

        for (i=0; i < nr; i++) {
            z_ = cr[i*2  ];
            x_ = cr[i*2+1];

            if (z_ < z[0] || z_ > z[Nz-1] ||
                x_ < x[0] || x_ > x[Nx-1])
                exit_error("receiver %d at x=%f z=%f outside domain",i,x_,z_);
        }
        
        printf("process %3d : receiver %4d at x=%8.2f z=%8.2f\n",mpi_rank,0,cr[1],cr[0]);
        printf("process %3d : receiver %4d at x=%8.2f z=%8.2f\n",mpi_rank,nr-1,cr[2*nr-1],cr[2*nr-2]);
    }

    /* transpose */
    if(backward) {
        data = (float *)calloc(nrt = nr * nt,sizeof(float)); // nr cutted
    
        for     (i=0; i < nr; i++)
            for (j=0; j < nt; j++)
                data[j*nr + i] = datat[i*nt + j];

        free(datat);
    } else if (savedata)
        data = (float *)calloc(nrt = nr * nt,sizeof(float));
}

void update(int back)
{
    int i,j,k,it,ic,ir,is,ih,ix,iz,iit,nrt_w,Ncz,Nrz,Ncx,Nrx,n[2],l[2],h[2],N[2];
    float dt2,o[2],d[2],*velo2,*alfa,*tmp,*pa,*po,*pb,*ua,*uo,*ub,*wa,*wo,*wb,*RZ,*RX;
    sf_complex tt,*CZ,*CX,*DZ,*DX;
    fftwf_plan fwdz,invz,fwdx,invx;
    FILE *Fuus,*Fimg=NULL,*Fcig=NULL,*Fdata;
#ifdef DEBUG
    FILE *Fuur;
#endif
    MPI_Status status;
    Int2 *I2;

    N[0] = Nz; d[0] = dz; o[0] = z[0];
    N[1] = Nx; d[1] = dx; o[1] = x[0];

    I2 = int2_init(N,o,d);

    dt2 = dt * 2.;

    n[0] = nz;  n[1] = nx;
    l[0] = bzl; l[1] = bxl;
    h[0] = bzh; h[1] = bxh;

    velo2= (float *)calloc(Nxz,sizeof(float));

    expand(velo2,velo,n,l,h);

    alfa = (float *)calloc(Nxz,sizeof(float));

    for (i=0; i < Nxz; i++) {
        alfa[i] = dt2 * velo2[i] * velo2[i];
    }

    Nrz = 2 * (Ncz = Nz/2 + 1);
    Nrx = 2 * (Ncx = Nx/2 + 1);

    RZ = (float *)fftwf_malloc(sizeof(float) * Nx * Nrz); CZ = (sf_complex *)RZ;
    RX = (float *)fftwf_malloc(sizeof(float) * Nz * Nrx); CX = (sf_complex *)RX;

    fwdz = fftwf_plan_dft_r2c_1d(Nz,RZ,(fftwf_complex *) CZ,FFTW_PATIENT);
    invz = fftwf_plan_dft_c2r_1d(Nz,(fftwf_complex *) CZ,RZ,FFTW_PATIENT);
    fwdx = fftwf_plan_dft_r2c_1d(Nx,RX,(fftwf_complex *) CX,FFTW_PATIENT);
    invx = fftwf_plan_dft_c2r_1d(Nx,(fftwf_complex *) CX,RX,FFTW_PATIENT);
    if (NULL == fwdz || NULL == invz ||
        NULL == fwdx || NULL == invx)
	exit_error("fftw planning failed");

    DZ = (sf_complex *)calloc(Ncz,sizeof(sf_complex));
    DX = (sf_complex *)calloc(Ncx,sizeof(sf_complex));
	
    for (tt = sf_cmplx(0.,2./Nz * SF_PI/dz), iz=0; iz < Ncz; iz++) DZ[iz] = kz[iz] * tt;
    for (tt = sf_cmplx(0.,2./Nx * SF_PI/dx), ix=0; ix < Ncx; ix++) DX[ix] = kx[ix] * tt;

    if (back) {
	if (NULL == (Fuus = fopen(fname_uus,"r")))
            exit_error("failed to open %s",fname_uus);
	if (NULL == (Fimg = fopen(fname_img,"w+")))
            exit_error("failed to open %s",fname_img);
	if (NULL == (Fcig = fopen(fname_cig,"w+")))
            exit_error("failed to open %s",fname_cig);
#ifdef DEBUG
	if (NULL == (Fuur = fopen(fname_uur,"w+")))
            exit_error("failed to open %s",fname_uur);
#endif

	uus = (float *)calloc(Nxz,sizeof(float));
	img = (float *)calloc(Nxz,sizeof(float));
	cig = (float *)calloc(NHCz,sizeof(float));

	fseek(Fuus,0,SEEK_END);

	for (i=0; i < Nxz;  i++) img[i] = 0.;
	for (i=0; i < NHCz; i++) cig[i] = 0.;
    } else {
	if (NULL == (Fuus = fopen(fname_uus,"w+")))
            exit_error("failed to open %s",fname_uus);
    }

    pa = (float *)calloc(Nxz,sizeof(float));
    po = (float *)calloc(Nxz,sizeof(float));
    pb = (float *)calloc(Nxz,sizeof(float));
    ua = (float *)calloc(Nxz,sizeof(float));
    uo = (float *)calloc(Nxz,sizeof(float));
    ub = (float *)calloc(Nxz,sizeof(float));
    wa = (float *)calloc(Nxz,sizeof(float));
    wo = (float *)calloc(Nxz,sizeof(float));
    wb = (float *)calloc(Nxz,sizeof(float));
    for (i=0; i < Nxz; i++)
        pa[i] = po[i] = pb[i] = \
	    ua[i] = uo[i] = ub[i] = \
	    wa[i] = wo[i] = wb[i] = 0.;

    for (it = 0; it < nt; it++) {
	iit = back ? nt-1-it : it;
	if (!(iit % 250))
	    printf("process %d : it = %d\n",mpi_rank,iit);

        if (!(iit % j3)) {
	    if (back) {
		fseek(Fuus,-Nxz * sizeof(float),SEEK_CUR);
		if (Nxz != fread(uus,sizeof(float),Nxz,Fuus))
		    exit_error("failed reading %s",fname_uus);
		fseek(Fuus,-Nxz * sizeof(float),SEEK_CUR);

		/* imaging condition */
#pragma omp parallel for schedule(dynamic) private(i)
		for (i=0; i < Nxz; i++)
		    img[i] += po[i] * uus[i];
		fseek(Fimg,0,SEEK_SET);
		if (Nxz != fwrite(img,sizeof(float),Nxz,Fimg))
		    exit_error("failed writing %s",fname_img);

		/* compute cig */
		for     (ih=0; ih < nH; ih++)
		    for (ic=0, ir=C0+ih*dH, is=C0-ih*dH; ic < nC; ic++, ir+=dC, is+=dC) {
                        ir = (ir >= nx) ? nx-1 : ir;
                        is = (is < 0)   ? 0    : is;
#pragma omp parallel for schedule(dynamic) private(iz)
			for (iz=0; iz < Nz; iz++) {
                            cig[(ih*nC + ic)*Nz + iz] += po[(bxl+ir)*Nz + iz] * uus[(bxl+is)*Nz + iz];
			}
                    }
		fseek(Fcig,0,SEEK_SET);
		if (NHCz != fwrite(cig,sizeof(float),NHCz,Fcig))
		    exit_error("failed writing %s",fname_cig);

#ifdef DEBUG
		if (Nxz != fwrite(po,sizeof(float),Nxz,Fuur))
		    exit_error("failed writing %s",fname_uur);
#endif
	    } else {
		if (Nxz != fwrite(po,sizeof(float),Nxz,Fuus))
		    exit_error("failed writing %s",fname_uus);
	    }
	}

        /* update particle momentum */
#pragma omp parallel for schedule(dynamic) private(iz,ix)
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++)
                RX[iz * Nrx + ix] = RZ[ix * Nrz + iz] = po[ix * Nz + iz];
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],(fftwf_complex *)&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,(fftwf_complex *)&CZ[ic],&RZ[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],(fftwf_complex *)&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,(fftwf_complex *)&CX[ic],&RX[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                ub[i] = ua[i] + dt2 * RX[k];
                wb[i] = wa[i] + dt2 * RZ[j];
	    }
	}

        /* update stress */
#pragma omp parallel for schedule(dynamic) private(iz,ix,i)
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz + iz;
                RX[iz * Nrx + ix] = uo[i];
                RZ[ix * Nrz + iz] = wo[i];
            }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],(fftwf_complex *)&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,(fftwf_complex *)&CZ[ic],&RZ[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],(fftwf_complex *)&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,(fftwf_complex *)&CX[ic],&RX[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                pb[i] = pa[i] + alfa[i] * RX[k] + alfa[i] * RZ[j];
	    }
	}

        sponge_apply(pb,n,sponge);
        sponge_apply(ub,n,sponge);
        sponge_apply(wb,n,sponge);

        if (!back) {
            int2_inject(&sour[it],pb,1,cs,I2);
        }

	/* read/write data */
	if (back) {
	    int2_inject(&data[iit*nr],pb,nr,cr,I2);
	}
	else if (savedata)
	    int2_apply(&data[iit*nr],pb,nr,cr,I2);

	tmp = pa; pa = po; po = pb; pb = tmp;
	tmp = ua; ua = uo; uo = ub; ub = tmp;
	tmp = wa; wa = wo; wo = wb; wb = tmp;
    } /* it */

    if (!back && savedata) {
        /* transpose */
        datat = (float *)calloc(nrt,sizeof(float));
        for     (ir=0; ir < nr; ir++)
            for (it=0; it < nt; it++)
                datat[ir*nt + it] = data[it*nr + ir];

        if (MASTER == mpi_rank) {
            if (NULL == (Fdata = fopen(fname_data,"w+")))
                exit_error("failed to open %s",fname_data);

            fwrite(datat,sizeof(float),nrt,Fdata);
            for (i=1; i < mpi_size; i++) {
		/* nrt_w <= nrt because of cutting */
		MPI_Recv(&nrt_w,1,MPI_INT,i,TAG_NR,MPI_COMM_WORLD,&status);

                MPI_Recv(datat,nrt_w,MPI_FLOAT,i,TAG_DA,MPI_COMM_WORLD,&status);
                MPI_Get_count(&status,MPI_FLOAT,&j);
                if (nrt_w != j) exit_error("failed receiving data. %d expected %d received",nrt_w,j);
		for (j=nrt_w; j < nrt; j++) datat[j] = 0.;
                fwrite(datat,sizeof(float),nrt,Fdata);
            }
        } else {
	    MPI_Send(&nrt,1,MPI_INT,MASTER,TAG_NR,MPI_COMM_WORLD);
            MPI_Send(datat,nrt,MPI_FLOAT,MASTER,TAG_DA,MPI_COMM_WORLD);
        }
	free(datat);
    }

    /* cleanup */
    fftwf_destroy_plan(fwdz);	fftwf_destroy_plan(invz);
    fftwf_destroy_plan(fwdx);	fftwf_destroy_plan(invx);
    fftwf_free(RZ); free(DZ);
    fftwf_free(RX); free(DX);
    free(pa); free(po); free(pb);
    free(ua); free(uo); free(ub);
    free(wa); free(wo); free(wb);
    free(velo2); free(alfa);
    if (fclose(Fuus)) exit_error("failed closing %s",fname_uus);
    if (back) {
	free(uus); free(img);
	if (fclose(Fimg)) exit_error("failed closing %s",fname_img);
	if (fclose(Fcig)) exit_error("failed closing %s",fname_cig);
#ifdef DEBUG
	if (fclose(Fuur)) exit_error("failed closing %s",fname_uur);
#endif
    }
}

void update_tau(int back)
{
    int i,j,k,it,ic,ir,is,ih,ix,iz,iit,nrt_w,Ncz,Nrz,Ncx,Nrx,n[2],l[2],h[2],N[2];
    float dt2,o[2],d[2],*velo2,*vmap2,*sigm2,*alfa,*beta,*gama,*tmp,*pa,*po,*pb,*ua,*uo,*ub,*wa,*wo,*wb,*RZ,*RX;
    sf_complex tt,*CZ,*CX,*DZ,*DX;
    fftwf_plan fwdz,invz,fwdx,invx;
    FILE *Fuus,*Fimg=NULL,*Fcig=NULL,*Fdata;
#ifdef DEBUG
    FILE *Fuur;
#endif
    MPI_Status status;
    Int2 *I2;

    N[0] = Nz; d[0] = dz; o[0] = z[0];
    N[1] = Nx; d[1] = dx; o[1] = x[0];

    I2 = int2_init(N,o,d);

    dt2 = dt * 2.;

    n[0] = nz;  n[1] = nx;
    l[0] = bzl; l[1] = bxl;
    h[0] = bzh; h[1] = bxh;

    velo2= (float *)calloc(Nxz,sizeof(float));
    vmap2= (float *)calloc(Nxz,sizeof(float));
    sigm2= (float *)calloc(Nxz,sizeof(float));

    expand(velo2,velo,n,l,h);
    expand(vmap2,vmap,n,l,h);
    expand(sigm2,sigm,n,l,h);

    alfa = (float *)calloc(Nxz,sizeof(float));
    beta = (float *)calloc(Nxz,sizeof(float));
    gama = (float *)calloc(Nxz,sizeof(float));

    for (i=0; i < Nxz; i++) {
        alfa[i] = dt2 * velo2[i] * velo2[i] / vmap2[i];
        beta[i] = dt2 * sigm2[i];
        gama[i] = dt2 * (1./vmap2[i]*1./vmap2[i] + sigm2[i]*sigm2[i]);
    }

    Nrz = 2 * (Ncz = Nz/2 + 1);
    Nrx = 2 * (Ncx = Nx/2 + 1);

    RZ = (float *)fftwf_malloc(sizeof(float) * Nx * Nrz); CZ = (sf_complex *)RZ;
    RX = (float *)fftwf_malloc(sizeof(float) * Nz * Nrx); CX = (sf_complex *)RX;

    fwdz = fftwf_plan_dft_r2c_1d(Nz,RZ,(fftwf_complex *)CZ,FFTW_PATIENT);
    invz = fftwf_plan_dft_c2r_1d(Nz,(fftwf_complex *)CZ,RZ,FFTW_PATIENT);
    fwdx = fftwf_plan_dft_r2c_1d(Nx,RX,(fftwf_complex *)CX,FFTW_PATIENT);
    invx = fftwf_plan_dft_c2r_1d(Nx,(fftwf_complex *)CX,RX,FFTW_PATIENT);
    if (NULL == fwdz || NULL == invz ||
        NULL == fwdx || NULL == invx)
	exit_error("fftw planning failed");

    DZ = (sf_complex *)calloc(Ncz,sizeof(sf_complex));
    DX = (sf_complex *)calloc(Ncx,sizeof(sf_complex));
	
    for (tt = sf_cmplx(0.,2./Nz * SF_PI/dz), iz=0; iz < Ncz; iz++) DZ[iz] = kz[iz] * tt;
    for (tt = sf_cmplx(0.,2./Nx * SF_PI/dx), ix=0; ix < Ncx; ix++) DX[ix] = kx[ix] * tt;

    if (back) {
	if (NULL == (Fuus = fopen(fname_uus,"r")))
            exit_error("failed to open %s",fname_uus);
	if (NULL == (Fimg = fopen(fname_img,"w+")))
            exit_error("failed to open %s",fname_img);
	if (NULL == (Fcig = fopen(fname_cig,"w+")))
            exit_error("failed to open %s",fname_cig);
#ifdef DEBUG
	if (NULL == (Fuur = fopen(fname_uur,"w+")))
            exit_error("failed to open %s",fname_uur);
#endif

	uus = (float *)calloc(Nxz,sizeof(float));
	img = (float *)calloc(Nxz,sizeof(float));
	cig = (float *)calloc(NHCz,sizeof(float));

	fseek(Fuus,0,SEEK_END);

	for (i=0; i < Nxz;  i++) img[i] = 0.;
	for (i=0; i < NHCz; i++) cig[i] = 0.;
    } else {
	if (NULL == (Fuus = fopen(fname_uus,"w+")))
            exit_error("failed to open %s",fname_uus);
    }

    pa = (float *)calloc(Nxz,sizeof(float));
    po = (float *)calloc(Nxz,sizeof(float));
    pb = (float *)calloc(Nxz,sizeof(float));
    ua = (float *)calloc(Nxz,sizeof(float));
    uo = (float *)calloc(Nxz,sizeof(float));
    ub = (float *)calloc(Nxz,sizeof(float));
    wa = (float *)calloc(Nxz,sizeof(float));
    wo = (float *)calloc(Nxz,sizeof(float));
    wb = (float *)calloc(Nxz,sizeof(float));
    for (i=0; i < Nxz; i++)
        pa[i] = po[i] = pb[i] = \
	    ua[i] = uo[i] = ub[i] = \
	    wa[i] = wo[i] = wb[i] = 0.;

    for (it = 0; it < nt; it++) {
	iit = back ? nt-1-it : it;
	if (!(iit % 250))
	    printf("process %d : it = %d\n",mpi_rank,iit);

        if (!(iit % j3)) {
	    if (back) {
		fseek(Fuus,-Nxz * sizeof(float),SEEK_CUR);
		if (Nxz != fread(uus,sizeof(float),Nxz,Fuus))
		    exit_error("failed reading %s",fname_uus);
		fseek(Fuus,-Nxz * sizeof(float),SEEK_CUR);

		/* imaging condition */
#pragma omp parallel for schedule(dynamic) private(i)
		for (i=0; i < Nxz; i++)
		    img[i] += po[i] * uus[i];
		fseek(Fimg,0,SEEK_SET);
		if (Nxz != fwrite(img,sizeof(float),Nxz,Fimg))
		    exit_error("failed writing %s",fname_img);

		/* compute cig */
		for     (ih=0; ih < nH; ih++)
		    for (ic=0, ir=C0+ih*dH, is=C0-ih*dH; ic < nC; ic++, ir+=dC, is+=dC) {
                        ir = (ir >= nx) ? nx-1 : ir;
                        is = (is < 0)   ? 0    : is;
#pragma omp parallel for schedule(dynamic) private(iz)
			for (iz=0; iz < Nz; iz++) {
                            cig[(ih*nC + ic)*Nz + iz] += po[(bxl+ir)*Nz + iz] * uus[(bxl+is)*Nz + iz];
			}
                    }
		fseek(Fcig,0,SEEK_SET);
		if (NHCz != fwrite(cig,sizeof(float),NHCz,Fcig))
		    exit_error("failed writing %s",fname_cig);

#ifdef DEBUG
		if (Nxz != fwrite(po,sizeof(float),Nxz,Fuur))
		    exit_error("failed writing %s",fname_uur);
#endif
	    } else {
		if (Nxz != fwrite(po,sizeof(float),Nxz,Fuus))
		    exit_error("failed writing %s",fname_uus);
	    }
	}

        /* update particle momentum */
#pragma omp parallel for schedule(dynamic) private(iz,ix)
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++)
                RX[iz * Nrx + ix] = RZ[ix * Nrz + iz] = po[ix * Nz + iz];
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],(fftwf_complex *) &CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,(fftwf_complex *) &CZ[ic],&RZ[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],(fftwf_complex *) &CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,(fftwf_complex *) &CX[ic],&RX[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                ub[i] = ua[i] +     dt2 * RX[k] + beta[i] * RZ[j];
                wb[i] = wa[i] + beta[i] * RX[k] + gama[i] * RZ[j];
	    }
	}

        /* update stress */
#pragma omp parallel for schedule(dynamic) private(iz,ix,i)
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz + iz;
                RX[iz * Nrx + ix] = uo[i] * vmap2[i];
                RZ[ix * Nrz + iz] = wo[i] * vmap2[i];
            }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],(fftwf_complex *) &CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,(fftwf_complex *) &CZ[ic],&RZ[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],(fftwf_complex *) &CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,(fftwf_complex *) &CX[ic],&RX[ir]);
        }
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                pb[i] = pa[i] + alfa[i] * RX[k] + alfa[i] * RZ[j];
	    }
	}

        sponge_apply(pb,n,sponge);
        sponge_apply(ub,n,sponge);
        sponge_apply(wb,n,sponge);

        if (!back) {
            int2_inject(&sour[it],pb,1,cs,I2);
        }

	/* read/write data */
	if (back) {
	    int2_inject(&data[iit*nr],pb,nr,cr,I2);
	}
	else if (savedata)
	    int2_apply(&data[iit*nr],pb,nr,cr,I2);

	tmp = pa; pa = po; po = pb; pb = tmp;
	tmp = ua; ua = uo; uo = ub; ub = tmp;
	tmp = wa; wa = wo; wo = wb; wb = tmp;
    } /* it */

    if (!back && savedata) {
        /* transpose */
        datat = (float *)calloc(nrt,sizeof(float));
        for     (ir=0; ir < nr; ir++)
            for (it=0; it < nt; it++)
                datat[ir*nt + it] = data[it*nr + ir];

        if (MASTER == mpi_rank) {
            if (NULL == (Fdata = fopen(fname_data,"w+")))
                exit_error("failed to open %s",fname_data);

            fwrite(datat,sizeof(float),nrt,Fdata);
            for (i=1; i < mpi_size; i++) {
		/* nrt_w <= nrt because of cutting */
		MPI_Recv(&nrt_w,1,MPI_INT,i,TAG_NR,MPI_COMM_WORLD,&status);

                MPI_Recv(datat,nrt_w,MPI_FLOAT,i,TAG_DA,MPI_COMM_WORLD,&status);
                MPI_Get_count(&status,MPI_FLOAT,&j);
                if (nrt_w != j) exit_error("failed receiving data. %d expected %d received",nrt_w,j);
		for (j=nrt_w; j < nrt; j++) datat[j] = 0.;
                fwrite(datat,sizeof(float),nrt,Fdata);
            }
        } else {
	    MPI_Send(&nrt,1,MPI_INT,MASTER,TAG_NR,MPI_COMM_WORLD);
            MPI_Send(datat,nrt,MPI_FLOAT,MASTER,TAG_DA,MPI_COMM_WORLD);
        }
	free(datat);
    }

    /* cleanup */
    fftwf_destroy_plan(fwdz);	fftwf_destroy_plan(invz);
    fftwf_destroy_plan(fwdx);	fftwf_destroy_plan(invx);
    fftwf_free(RZ); free(DZ);
    fftwf_free(RX); free(DX);
    free(pa); free(po); free(pb);
    free(ua); free(uo); free(ub);
    free(wa); free(wo); free(wb);
    free(velo2); free(alfa);
    if (fclose(Fuus)) exit_error("failed closing %s",fname_uus);
    if (back) {
	free(uus); free(img);
	if (fclose(Fimg)) exit_error("failed closing %s",fname_img);
	if (fclose(Fcig)) exit_error("failed closing %s",fname_cig);
#ifdef DEBUG
	if (fclose(Fuur)) exit_error("failed closing %s",fname_uur);
#endif
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    init(argc,argv);

    if (forward)  {
        if (tau) update_tau(0);
        else     update(0);
    }
    if (backward) {
        if (tau) update_tau(1);
        else     update(1);
    }

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
