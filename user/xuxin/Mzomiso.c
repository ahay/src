/* zero-offset isotropic reverse-time migration
 * exploding reflector modeling : < imag.rsf sfzomiso inv=y > data.rsf
 * zero-offset migration        : < data.rsf sfzomiso inv=n > imag.rsf
 * forward modeling             : < data.rsf sfzomiso inv=n > imag.rsf 
 Need (1) velo *= 2 (2) nr=1 (3) imag.rsf is useless */

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

/* -------- Equations (cartesian) ------------- */
/*                                              */
/* dp/dt = v^2 * (du/dx + dw/dz)                */
/* du/dt = dp/dx                                */
/* dw/dt = dp/dz                                */
/*                                              */
/* -------- Equations (nonortho tau) ---------- */
/*                                              */
/* dp/dt = v^2/vm * d/dx (vm * u)               */
/*       + v^2/vm * d/dz (vm * w)               */
/* du/dt =      dp/dx +              si * dp/dz */
/* dw/dt = si * dp/dx + (si^2 + 1/vm^2) * dp/dz */
/*                                              */
/* p     : stress                               */
/* u,w   : particle momentum * (-1)             */
/* v     : P-wave velocity                      */
/* vm    : mapping velocity                     */
/* si    : sigma                                */
/* z     : vertical time if tau=y               */
/* -------------------------------------------- */

/* -------------------------------------- */
/* File   Size          Description       */
/* -------------------------------------- */
/* velo=      [nx][nz]  velocity          */
/* cr=        [nr][2 ]  receiver position */
/* in/out=    [nx][nz]  reflectivity      */
/* out/in=    [nt][nr]  zero-offset data  */
/* wave=  [n3][Nx][Nz]  wavefield         */
/* vmap=      [nx][nz]  mapping velocity  */
/* sigm=      [nx][nz]  sigma             */
/* -------------------------------------- */

/* stability condition : dx / (dt*vmax/2) >= 5 */
/* avoid having 99% null traces in data, it causes high k artefacts. */
/* orthogonal tau (lateral constant) is not implemented */

#include <rsf.h>
#include <fftw3.h>
#include "util.h"
#include "sponge.h"
#include "interp.h"

#ifdef _OPENMP
#include <omp.h>
int omp_nth;
#endif

#define I _Complex_I

int nz; float dz,z0,*z,*kz; int bzl,bzh,Nz;
int nx; float dx,x0,*x,*kx; int bxl,bxh,Nx;
int nt; float dt,   *t;
int n3; float d3;           int j3;

int nr,nxz,Nxz,nrt;

Sponge *sponge;

sf_file Fvelo; float *velo; /*     [nx][nz] */
sf_file Fimag; float *imag; /*     [nx][nz] */
sf_file Fdata; float *data; /*     [nt][nr] */
sf_file Fcr;   float *cr;   /*     [nr][2 ] */
sf_file Fwave;              /* [n3][Nx][Nz] */

bool inv; /* true : modeling
             false: migration */

/* ============ tau parameters ============ */

bool tau; /* true : tau
             false: cartesian */
sf_file Fvmap; float *vmap; /*     [nx][nz] */
sf_file Fsigm; float *sigm; /*     [nx][nz] */

void init()
{
    int i,b[4];
    float c[4],czl,czh,cxl,cxh,z_,x_,tmpf,eps;
    bool opt,verb;

    if (!sf_getbool("inv",&inv))   inv = false; /* if y, modeling; if n, migration */
    if (!sf_getbool("tau",&tau))   tau = false; /* if y, tau domain; if n, cartesian */
    if (!sf_getbool("verb",&verb)) verb= false; /* verbosity */
    if (!sf_getbool("opt",&opt))   opt = false; /* optimze fft size */

    Fvelo = sf_input("velo");
    Fcr   = sf_input("cr");
    Fwave = sf_output("wave");
    if (inv) { Fimag = sf_input("in");   Fdata = sf_output("out"); }
    else     { Fimag = sf_output("out"); Fdata = sf_input("in");   }
    if (tau) { Fvmap = sf_input("vmap"); Fsigm = sf_input("sigm"); }

    if (!sf_histint  (Fvelo,"n1",&nz) ||
	!sf_histint  (Fvelo,"n2",&nx) ||
        !sf_histfloat(Fvelo,"d1",&dz) ||
        !sf_histfloat(Fvelo,"d2",&dx)) sf_error("Need n1= n2= d1= d2= in input");
    if (!sf_histfloat(Fvelo,"o1",&z0)) z0 = 0.;
    if (!sf_histfloat(Fvelo,"o2",&x0)) x0 = 0.;

    if (!sf_histint(Fcr,"n1",&i) || i != 2) sf_error("Need n1=2 in cr");
    if (!sf_histint(Fcr,"n2",&nr)) nr = 1;

    if (inv) {
        if (!sf_getint  ("nt",&nt)) nt = 1;  /* time n (if inv=y) */
        if (!sf_getfloat("dt",&dt)) dt = 1.; /* time d (if inv=y) */

        if (!sf_histint(Fimag,"n1",&i) || i != nz ||
            !sf_histint(Fimag,"n2",&i) || i != nx)
            sf_error("Need n1=%d n2=%d in input",nz,nx);

        sf_putint(Fdata,"n1",nr); sf_putint  (Fdata,"d1",1);  sf_putint  (Fdata,"o1",0);
        sf_putint(Fdata,"n2",nt); sf_putfloat(Fdata,"d2",dt); sf_putfloat(Fdata,"o2",0.);
        sf_putstring(Fdata,"label1","r"); sf_putstring(Fdata,"unit1","");
        sf_putstring(Fdata,"label2","t"); sf_putstring(Fdata,"unit2","sec");
    } else {
        if (!sf_histint  (Fdata,"n1",&i ) || i != nr ||
            !sf_histint  (Fdata,"n2",&nt) ||
            !sf_histfloat(Fdata,"d2",&dt))
            sf_error("Need n1=%d n2= d2= in input",nr);
        
        sf_putint(Fimag,"n1",nz); sf_putfloat(Fimag,"d1",dz); sf_putfloat(Fimag,"o1",z0);
        sf_putint(Fimag,"n2",nx); sf_putfloat(Fimag,"d2",dx); sf_putfloat(Fimag,"o2",x0);
        sf_putstring(Fimag,"label1",tau ? "tau" : "z");
        sf_putstring(Fimag,"unit1",tau ? "sec" : "m");
        sf_putstring(Fimag,"label2","x"); sf_putstring(Fimag,"unit2","m");
    }

    if (!sf_getint  ("bzl",&bzl)) bzl = 0;
    if (!sf_getint  ("bzh",&bzh)) bzh = 0;
    if (!sf_getint  ("bxl",&bxl)) bxl = 0;
    if (!sf_getint  ("bxh",&bxh)) bxh = 0;
    if (!sf_getfloat("czl",&czl)) czl = 1.;
    if (!sf_getfloat("czh",&czh)) czh = 1.;
    if (!sf_getfloat("cxl",&cxl)) cxl = 1.;
    if (!sf_getfloat("cxh",&cxh)) cxh = 1.;
    if (!sf_getint  ("n3",&n3))   n3 = nt; /* wave time n */
    if (tau && !sf_getfloat("eps",&eps)) eps = 1; /* regularize sigma */

    if ((nt-1) % (n3-1)) sf_warning("j3 not round nt=%d n3=%d",nt,n3);
    j3 = (nt-1) / (n3-1);
    d3 = (inv ? 1. : -1.) * dt * j3;

    if (opt) {
	fft_expand(nz,&bzl,&bzh);
	fft_expand(nx,&bxl,&bxh);
    }

    Nz = nz + bzl + bzh;
    Nx = nx + bxl + bxh;

    z = sf_floatalloc(Nz);
    x = sf_floatalloc(Nx);
    t = sf_floatalloc(nt);

    for (i=0; i < Nz; i++) z[i] = z0 + dz * (i-bzl);
    for (i=0; i < Nx; i++) x[i] = x0 + dx * (i-bxl);
    for (i=0; i < nt; i++) t[inv ? i : nt-1-i] = dt * i;

    sf_putint(Fwave,"n1",Nz); sf_putfloat(Fwave,"d1",dz); sf_putfloat(Fwave,"o1",z[0]);
    sf_putint(Fwave,"n2",Nx); sf_putfloat(Fwave,"d2",dx); sf_putfloat(Fwave,"o2",x[0]);
    sf_putint(Fwave,"n3",n3); sf_putfloat(Fwave,"d3",d3); sf_putfloat(Fwave,"o3",t[0]);
    sf_putstring(Fwave,"label1",tau ? "tau" : "z");
    sf_putstring(Fwave,"unit1",tau ? "sec" : "m");
    sf_putstring(Fwave,"label2","x"); sf_putstring(Fwave,"unit2","m");
    sf_putstring(Fwave,"label3","t"); sf_putstring(Fwave,"unit3","sec");

    nxz = nx * nz;
    Nxz = Nx * Nz;
    nrt = nr * nt;

    kz = sf_floatalloc(Nz); compute_k(kz,Nz);
    kx = sf_floatalloc(Nx);	compute_k(kx,Nx);

    b[0] = bzl; b[1] = bzh; b[2] = bxl; b[3] = bxh;
    c[0] = czl; c[1] = czh; c[2] = cxl; c[3] = cxh;
    sponge = sponge_init(b,c);

    velo = sf_floatalloc(nxz);
    imag = sf_floatalloc(nxz);
    data = sf_floatalloc(nrt);
    cr   = sf_floatalloc(nr*2);

    sf_floatread(velo,nxz,Fvelo);

    if (inv) sf_floatread(imag,nxz,Fimag);
    else     sf_floatread(data,nrt,Fdata);

    for (i=0; i < nr; i++) {
        sf_floatread(&cr[i*2],2,Fcr);

        z_ = cr[i*2  ];
        x_ = cr[i*2+1];

	if (z_ < z[0] || z_ > z[Nz-1] ||
            x_ < x[0] || x_ > x[Nx-1])
	    sf_error("receiver %d at x=%g z=%g outside domain",i,x_,z_);
    }

    if (tau) {
        vmap = sf_floatalloc(nxz);
        sigm = sf_floatalloc(nxz);
       
        if (!sf_histint  (Fvmap,"n1",&i) || i != nz ||
            !sf_histint  (Fvmap,"n2",&i) || i != nx)
            sf_error("Need n1=%d n2=%d in vmap",nz,nx);

        if (!sf_histint(Fsigm,"n1",&i) || i != nz ||
            !sf_histint(Fsigm,"n2",&i) || i != nx)
            sf_error("Need n1=%d n2=%d in sigm",nz,nx);

        sf_floatread(vmap,nxz,Fvmap);
        sf_floatread(sigm,nxz,Fsigm);

	/* regularize sigma */
	/* usually sigma at 1e-4 magnitude */
        for (i=0; i < nxz; i++)
            if (fabs(tmpf = sigm[i]) >= eps) sigm[i] = eps * (tmpf > 0 ? 1. : -1.);
    }

    if (verb) {
	sf_warning("%s domain isotropic %s",tau ? "tau" : "cartesian",inv ? "modeling" : "migration");
	sf_warning("receiver ir=%4d at x=%g z=%g",0,cr[1],cr[0]);
	sf_warning("receiver ir=%4d at x=%g z=%g",nr-1,cr[2*nr-1],cr[2*nr-2]);
	sf_warning("bzl=%d wzl[%d]=%g wzl[%d]=%g",bzl,0,sponge->wzl[0],sponge->bzl-1,sponge->wzl[bzl-1]);
	sf_warning("bzh=%d wzh[%d]=%g wzh[%d]=%g",bzh,0,sponge->wzh[0],sponge->bzh-1,sponge->wzh[bzh-1]);
	sf_warning("bxl=%d wxl[%d]=%g wxl[%d]=%g",bxl,0,sponge->wxl[0],sponge->bxl-1,sponge->wxl[bxl-1]);
	sf_warning("bxh=%d wxh[%d]=%g wxh[%d]=%g",bxh,0,sponge->wxh[0],sponge->bxh-1,sponge->wxh[bxh-1]);
	sf_warning("n3=%d j3=%d d3=%g",n3,j3,d3);
    }

#ifdef _OPENMP
#pragma omp parallel
    omp_nth = omp_get_num_threads();
    sf_warning("omp using %d threads",omp_nth);
#endif

    sf_fileclose(Fvelo);
    sf_fileclose(Fcr);
    if (tau) {sf_fileclose(Fvmap); sf_fileclose(Fsigm);}
}

void update()
{
    int i,j,k,it,iit,ic,ir,ix,iz,Ncz,Nrz,Ncx,Nrx,n[2],l[2],h[2],N[2];
    float dt2,o[2],d[2],*velo2,*alfa,*pa,*po,*pb,*ua,*uo,*ub,*wa,*wo,*wb,*tmp,*RZ,*RX;
    float complex tt,*CZ,*CX,*DZ,*DX;
    fftwf_plan fwdz,invz,fwdx,invx;
    Int2 *I2;

    N[0] = Nz; d[0] = dz; o[0] = z[0];
    N[1] = Nx; d[1] = dx; o[1] = x[0];
    
    I2 = int2_init(N,o,d);

    dt2 = dt * 2.;

    n[0] = nz;  n[1] = nx;
    l[0] = bzl; l[1] = bxl;
    h[0] = bzh; h[1] = bxh;

    velo2= sf_floatalloc(Nxz);

    expand(velo2,velo,n,l,h);


    for (i=0; i < Nxz; i++) velo2[i] /= 2.; /* half velocity */
    
    alfa = sf_floatalloc(Nxz);

    for (i=0; i < Nxz; i++) {
        alfa[i] = dt2 * velo2[i]*velo2[i];
    }

    Nrz = 2 * (Ncz = Nz/2 + 1);
    Nrx = 2 * (Ncx = Nx/2 + 1);

    RZ = (float *)fftwf_malloc(sizeof(float) * Nx * Nrz); CZ = (float complex *)RZ;
    RX = (float *)fftwf_malloc(sizeof(float) * Nz * Nrx); CX = (float complex *)RX;

    fwdz = fftwf_plan_dft_r2c_1d(Nz,RZ,CZ,FFTW_ESTIMATE);
    invz = fftwf_plan_dft_c2r_1d(Nz,CZ,RZ,FFTW_ESTIMATE);
    fwdx = fftwf_plan_dft_r2c_1d(Nx,RX,CX,FFTW_ESTIMATE);
    invx = fftwf_plan_dft_c2r_1d(Nx,CX,RX,FFTW_ESTIMATE);
    if (NULL == fwdz || NULL == invz ||
        NULL == fwdx || NULL == invx)
	sf_error("fftw planning failed");

    DZ = sf_complexalloc(Ncz);
    DX = sf_complexalloc(Ncx);

    for (tt=sf_cmplx(0.,2./Nz * SF_PI/dz), iz=0; iz < Ncz; iz++) DZ[iz] = kz[iz] * tt;
    for (tt=sf_cmplx(0.,2./Nx * SF_PI/dx), ix=0; ix < Ncx; ix++) DX[ix] = kx[ix] * tt;

    pa = sf_floatalloc(Nxz); po = sf_floatalloc(Nxz); pb = sf_floatalloc(Nxz);
    ua = sf_floatalloc(Nxz); uo = sf_floatalloc(Nxz); ub = sf_floatalloc(Nxz);
    wa = sf_floatalloc(Nxz); wo = sf_floatalloc(Nxz); wb = sf_floatalloc(Nxz);
    for (i=0; i < Nxz; i++) {
        pa[i] = po[i] = pb[i] = \
	    wa[i] = wo[i] = wb[i] = \
	    ua[i] = uo[i] = ub[i] = 0.;
    }
    
    /* load image */
    if (inv) {
        for    (ix=0; ix < nx; ix++)
            for(iz=0; iz < nz; iz++)
                pa[(ix+bxl) * Nz + iz+bzl] = \
		    po[(ix+bxl) * Nz + iz+bzl] = imag[ix * nz + iz];
    }

    for (it = 0; it < nt; it++) {
	sf_warning("it=%d;",iit = inv ? it : nt-1-it);

        /* update particle momentum */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix)
#endif 
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++)
                RX[iz * Nrx + ix] = RZ[ix * Nrz + iz] = po[ix * Nz + iz];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif 
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,&CZ[ic],&RZ[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif        
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,&CX[ic],&RX[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i)
#endif
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++) {
                i = ix * Nz + iz;
                RX[iz * Nrx + ix] = uo[i];
                RZ[ix * Nrz + iz] = wo[i];
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,&CZ[ic],&RZ[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,&CX[ic],&RX[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
#endif
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                pb[i] = pa[i] + alfa[i] * RX[k] + alfa[i] * RZ[j];
	    }
	}

        if (inv) {
            int2_apply(&data[iit*nr],pb,nr,cr,I2);
        } else {
            int2_inject(&data[iit*nr],pb,nr,cr,I2);
        }

        sponge_apply(pb,n,sponge);
        sponge_apply(ub,n,sponge);
        sponge_apply(wb,n,sponge);

	if (!(it % j3))
	    sf_floatwrite(pb,Nxz,Fwave);

	tmp = pa; pa = po; po = pb; pb = tmp;
	tmp = ua; ua = uo; uo = ub; ub = tmp;
	tmp = wa; wa = wo; wo = wb; wb = tmp;
    } /* it */
    sf_warning("\n");

    if (inv) {
        sf_floatwrite(data,nrt,Fdata);
    } else {
        for    (ix=0; ix < nx; ix++)
            for(iz=0; iz < nz; iz++)
                imag[ix * nz + iz] = po[(ix+bxl) * Nz + iz+bzl];

        sf_floatwrite(imag,nxz,Fimag);
    }

    sf_fileclose(Fwave);

    free(pa);    free(po);    free(pb);
    free(ua);    free(uo);    free(ub);
    free(wa);    free(wo);    free(wb);
    free(velo2);
    free(alfa);
    fftwf_free(RZ); free(DZ);
    fftwf_free(RX); free(DX);
    fftwf_destroy_plan(fwdz); fftwf_destroy_plan(invz);
    fftwf_destroy_plan(fwdx); fftwf_destroy_plan(invx);
    int2_free(I2);
}

void update_tau()
{
    int i,j,k,it,iit,ic,ir,ix,iz,Ncz,Nrz,Ncx,Nrx,n[2],l[2],h[2],N[2];
    float dt2,o[2],d[2],*velo2,*vmap2,*sigm2,*alfa,*beta,*gama,*pa,*po,*pb,*ua,*uo,*ub,*wa,*wo,*wb,*tmp,*RZ,*RX;
    float complex tt,*CZ,*CX,*DZ,*DX;
    fftwf_plan fwdz,invz,fwdx,invx;
    Int2 *I2;

    N[0] = Nz; d[0] = dz; o[0] = z[0];
    N[1] = Nx; d[1] = dx; o[1] = x[0];

    I2 = int2_init(N,o,d);

    dt2 = dt * 2.;

    n[0] = nz;  n[1] = nx;
    l[0] = bzl; l[1] = bxl;
    h[0] = bzh; h[1] = bxh;

    velo2= sf_floatalloc(Nxz);
    vmap2= sf_floatalloc(Nxz);
    sigm2= sf_floatalloc(Nxz);

    expand(velo2,velo,n,l,h);
    expand(vmap2,vmap,n,l,h);
    expand(sigm2,sigm,n,l,h);

    for (i=0; i < Nxz; i++) velo2[i] /= 2.;

    alfa = sf_floatalloc(Nxz);
    beta = sf_floatalloc(Nxz);
    gama = sf_floatalloc(Nxz);

    for (i=0; i < Nxz; i++) {
        alfa[i] = dt2 * velo2[i] * velo2[i] / vmap2[i];
        beta[i] = dt2 * sigm2[i];
        gama[i] = dt2 * (1./vmap2[i] * 1./vmap2[i] + sigm2[i] * sigm2[i]);
    }

    Nrz = 2 * (Ncz = Nz/2 + 1);
    Nrx = 2 * (Ncx = Nx/2 + 1);

    RZ = (float *)fftwf_malloc(sizeof(float) * Nx * Nrz); CZ = (float complex *)RZ;
    RX = (float *)fftwf_malloc(sizeof(float) * Nz * Nrx); CX = (float complex *)RX;

    fwdz = fftwf_plan_dft_r2c_1d(Nz,RZ,CZ,FFTW_ESTIMATE);
    invz = fftwf_plan_dft_c2r_1d(Nz,CZ,RZ,FFTW_ESTIMATE);
    fwdx = fftwf_plan_dft_r2c_1d(Nx,RX,CX,FFTW_ESTIMATE);
    invx = fftwf_plan_dft_c2r_1d(Nx,CX,RX,FFTW_ESTIMATE);
    if (NULL == fwdz || NULL == invz ||
        NULL == fwdx || NULL == invx)
	sf_error("fftw planning failed");

    DZ = sf_complexalloc(Ncz);
    DX = sf_complexalloc(Ncx);

    for (tt=sf_cmplx(0.,2./Nz * SF_PI/dz), iz=0; iz < Ncz; iz++) DZ[iz] = kz[iz] * tt;
    for (tt=sf_cmplx(0.,2./Nx * SF_PI/dx), ix=0; ix < Ncx; ix++) DX[ix] = kx[ix] * tt;

    pa = sf_floatalloc(Nxz); po = sf_floatalloc(Nxz); pb = sf_floatalloc(Nxz);
    ua = sf_floatalloc(Nxz); uo = sf_floatalloc(Nxz); ub = sf_floatalloc(Nxz);
    wa = sf_floatalloc(Nxz); wo = sf_floatalloc(Nxz); wb = sf_floatalloc(Nxz);
    for (i=0; i < Nxz; i++) {
        pa[i] = po[i] = pb[i] = \
	    wa[i] = wo[i] = wb[i] = \
	    ua[i] = uo[i] = ub[i] = 0.;
    }
    
    /* load image */
    if (inv) {
        for    (ix=0; ix < nx; ix++)
            for(iz=0; iz < nz; iz++)
                pa[(ix+bxl) * Nz + iz+bzl] = \
		    po[(ix+bxl) * Nz + iz+bzl] = imag[ix * nz + iz];
    }

    for (it = 0; it < nt; it++) {
	sf_warning("it=%d;",iit = inv ? it : nt-1-it);

        /* update particle momentum */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix)
#endif
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++)
                RX[iz * Nrx + ix] = RZ[ix * Nrz + iz] = po[ix * Nz + iz];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,&CZ[ic],&RZ[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif        
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,&CX[ic],&RX[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i)
#endif
	for     (ix=0; ix < Nx; ix++)
	    for (iz=0; iz < Nz; iz++) {
                i = ix * Nz + iz;
                RX[iz * Nrx + ix] = uo[i] * vmap2[i];
                RZ[ix * Nrz + iz] = wo[i] * vmap2[i];
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif
        for (ix=0; ix < Nx; ix++) {
            ic = ix * Ncz;
            ir = ix * Nrz;
            fftwf_execute_dft_r2c(fwdz,&RZ[ir],&CZ[ic]);
            for (iz=0; iz < Ncz; iz++)
                CZ[ix*Ncz + iz] *= DZ[iz];
            fftwf_execute_dft_c2r(invz,&CZ[ic],&RZ[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,ic,ir)
#endif
        for (iz=0; iz < Nz; iz++) {
            ic = iz * Ncx;
            ir = iz * Nrx;
            fftwf_execute_dft_r2c(fwdx,&RX[ir],&CX[ic]);
            for (ix=0; ix < Ncx; ix++)
                CX[iz*Ncx + ix] *= DX[ix];
            fftwf_execute_dft_c2r(invx,&CX[ic],&RX[ir]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ix,i,j,k)
#endif
	for     (ix=0; ix < Nx; ix++) {
	    for (iz=0; iz < Nz; iz++) {
		i = ix * Nz  + iz;
		j = ix * Nrz + iz;
                k = iz * Nrx + ix;

                pb[i] = pa[i] + alfa[i] * RX[k] + alfa[i] * RZ[j];
	    }
	}

        if (inv) {
            int2_apply(&data[iit*nr],pb,nr,cr,I2);
        } else {
            int2_inject(&data[iit*nr],pb,nr,cr,I2);
        }

        sponge_apply(pb,n,sponge);
        sponge_apply(ub,n,sponge);
        sponge_apply(wb,n,sponge);

	if (!(it % j3))
	    sf_floatwrite(pb,Nxz,Fwave);

	tmp = pa; pa = po; po = pb; pb = tmp;
	tmp = ua; ua = uo; uo = ub; ub = tmp;
	tmp = wa; wa = wo; wo = wb; wb = tmp;
    } /* it */
    sf_warning("\n");

    if (inv) {
        sf_floatwrite(data,nrt,Fdata);
    } else {
        for    (ix=0; ix < nx; ix++)
            for(iz=0; iz < nz; iz++)
                imag[ix * nz + iz] = po[(ix+bxl) * Nz + iz+bzl];

        sf_floatwrite(imag,nxz,Fimag);
    }

    sf_fileclose(Fwave);

    free(pa);    free(po);    free(pb);
    free(ua);    free(uo);    free(ub);
    free(wa);    free(wo);    free(wb);
    free(velo2); free(vmap2); free(sigm2);
    free(alfa);  free(beta);  free(gama);
    fftwf_free(RZ); free(DZ);
    fftwf_free(RX); free(DX);
    fftwf_destroy_plan(fwdz); fftwf_destroy_plan(invz);
    fftwf_destroy_plan(fwdx); fftwf_destroy_plan(invx);
    int2_free(I2);
}

int main(int argc, char *argv[])
{
    sf_init(argc,argv);

    init();

    if (tau) update_tau();
    else     update();

    return 0;
}
