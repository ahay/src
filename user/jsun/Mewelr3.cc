/* Lowrank symbol approxiamtion for 3-D recursive integral time extrapolation of elastic waves */
/*
   Copyright (C) 2016 The University of Texas at Austin 
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

/* standard header */
#include <rsf.hh>
#include <rsf.h>
#include <time.h>
#include <assert.h>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

/* automatically generated header files */
extern "C"{
#include "eigen3x3.h"
}

/* static variables */
static float dt;
static bool tric,tstp,pseu;
static std::valarray<float>  C11,C12,C13,C22,C23,C33,C44,C55,C66,Q1,Q2;
static std::valarray<float>  C14,C15,C16,C24,C25,C26,C34,C35,C36,C45,C46,C56;
static std::valarray<double> rkz, rkx, rky;
static int component, mode;

/* subroutines */
static int sample(vector<int>& rs, vector<int>& cs, ZpxNumMat& res);

void expand3d(std::valarray<float>  a, std::valarray<float>& b, int nb, int nz, int nx, int ny);

/* global variables needed by lapack SVD routines */
char    jobz='V';
char    uplo='U';
int     M=3;
int     LDA=M;
int     LWORK=8*M;
int     INFO;
double  Chr[9], ww[3], work[24];

/*****************************************************************************************/
int main(int argc, char* argv[])
{
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* execution parameters                                       */
    /*------------------------------------------------------------*/

    /* wall clock */
    clock_t t1, t2, t3;
    float   timespent;
    t1=clock();

    iRSF par(0);
    bool verb;
    par.get("verb",verb,false);      // verbosity flag
    bool tilt;
    par.get("tilt",tilt,false);      // tilting of TTI
    int seed;
    par.get("tric",tric,false);      // triclinic anisotropy
    par.get("tstp",tstp,false);      // twostep propagator
    par.get("pseu",pseu,false);      // pseudo-spectral propagator
    par.get("seed",seed,time(NULL)); // seed for random number generator
    //srand48(seed);
    par.get("mode",mode,0);          // mode of decomposition: 0->mixed, 1->p, 2->s
    int jump;
    par.get("jump",jump,1);          // jump step for reduced lowrank decomposition
    float eps;
    par.get("eps",eps,1.e-6);        // tolerance
    int npk;
    par.get("npk",npk,20);           // maximum rank
    int nb;
    par.get("nb",nb,0);              // boundary padding
    par.get("dt",dt,1.e-3);          // time step size

    /*------------------------------------------------------------*/
    /* setup I/O files                                            */
    /*------------------------------------------------------------*/

    /* orthorhombic stiffness tensor (Voigt recipe) */
    /* c11 c12 c13 .   .   . 
       .   c22 c23 .   .   . 
       .   .   c33 .   .   . 
       .   .   .   c44 .   .
       .   .   .   .   c55 .
       .   .   .   .   .   c66 */

    iRSF c11,c12("c12"),c13("c13"),c22("c22"),c23("c23"),c33("c33"),c44("c44"),c55("c55"),c66("c66");

    /* read/write axes */
    int   nx, ny, nz;
    float dx, dy, dz;
    float ox, oy, oz;
    c11.get("n1",nz); c11.get("d1",dz); c11.get("o1",oz);
    c11.get("n2",nx); c11.get("d2",dx); c11.get("o2",ox);
    c11.get("n3",ny); c11.get("d3",dy); c11.get("o3",oy);

    /* calculate padded dimension */
    int nzpad = nz+2*nb;
    int nxpad = nx+2*nb;
    int nypad = ny+2*nb;
    float ozpad = oz-nb*dz;
    float oxpad = ox-nb*dx;
    float oypad = oy-nb*dy;
    int nxyz = nzpad*nxpad*nypad;

    /* read in parameter arrays */
    if (nb==0) {
        C11.resize(nxyz); c11>>C11;
        C12.resize(nxyz); c12>>C12;
        C13.resize(nxyz); c13>>C13;
        C22.resize(nxyz); c22>>C22;
        C23.resize(nxyz); c23>>C23;
        C33.resize(nxyz); c33>>C33;
        C44.resize(nxyz); c44>>C44;
        C55.resize(nxyz); c55>>C55;
        C66.resize(nxyz); c66>>C66;
        if (tric) {
            iRSF c14("c14"), c15("c15"), c16("c16"), 
                 c24("c24"), c25("c25"), c26("c26"), 
                 c34("c34"), c35("c35"), c36("c36"), 
                 c45("c45"), c46("c46"), c56("c56"); 
            C14.resize(nxyz); c14>>C14;
            C15.resize(nxyz); c15>>C15;
            C16.resize(nxyz); c16>>C16;
            C24.resize(nxyz); c24>>C24;
            C25.resize(nxyz); c25>>C25;
            C26.resize(nxyz); c26>>C26;
            C34.resize(nxyz); c34>>C34;
            C35.resize(nxyz); c35>>C35;
            C36.resize(nxyz); c36>>C36;
            C45.resize(nxyz); c45>>C45;
            C46.resize(nxyz); c46>>C46;
            C56.resize(nxyz); c56>>C56;
        }
        if (tilt) {
            iRSF q1("theta"), q2("phi");
            Q1.resize(nxyz); q1>>Q1;
            Q2.resize(nxyz); q2>>Q2;
        } else {
            Q1.resize(nxyz,0.);
            Q2.resize(nxyz,0.);
        }
    } else {
        std::valarray<float> tt(nz*nx*ny);
        C11.resize(nxyz); c11>>tt; expand3d(tt,C11,nb,nz,nx,ny);
        C12.resize(nxyz); c12>>tt; expand3d(tt,C12,nb,nz,nx,ny);
        C13.resize(nxyz); c13>>tt; expand3d(tt,C13,nb,nz,nx,ny);
        C22.resize(nxyz); c22>>tt; expand3d(tt,C22,nb,nz,nx,ny);
        C23.resize(nxyz); c23>>tt; expand3d(tt,C23,nb,nz,nx,ny);
        C33.resize(nxyz); c33>>tt; expand3d(tt,C33,nb,nz,nx,ny);
        C44.resize(nxyz); c44>>tt; expand3d(tt,C44,nb,nz,nx,ny);
        C55.resize(nxyz); c55>>tt; expand3d(tt,C55,nb,nz,nx,ny);
        C66.resize(nxyz); c66>>tt; expand3d(tt,C66,nb,nz,nx,ny);
        if (tric) {
            iRSF c14("c14"), c15("c15"), c16("c16"), 
                 c24("c24"), c25("c25"), c26("c26"), 
                 c34("c34"), c35("c35"), c36("c36"), 
                 c45("c45"), c46("c46"), c56("c56"); 
            C14.resize(nxyz); c14>>tt; expand3d(tt,C14,nb,nz,nx,ny);
            C15.resize(nxyz); c15>>tt; expand3d(tt,C15,nb,nz,nx,ny);
            C16.resize(nxyz); c16>>tt; expand3d(tt,C16,nb,nz,nx,ny);
            C24.resize(nxyz); c24>>tt; expand3d(tt,C24,nb,nz,nx,ny);
            C25.resize(nxyz); c25>>tt; expand3d(tt,C25,nb,nz,nx,ny);
            C26.resize(nxyz); c26>>tt; expand3d(tt,C26,nb,nz,nx,ny);
            C34.resize(nxyz); c34>>tt; expand3d(tt,C34,nb,nz,nx,ny);
            C35.resize(nxyz); c35>>tt; expand3d(tt,C35,nb,nz,nx,ny);
            C36.resize(nxyz); c36>>tt; expand3d(tt,C36,nb,nz,nx,ny);
            C45.resize(nxyz); c45>>tt; expand3d(tt,C45,nb,nz,nx,ny);
            C46.resize(nxyz); c46>>tt; expand3d(tt,C46,nb,nz,nx,ny);
            C56.resize(nxyz); c56>>tt; expand3d(tt,C56,nb,nz,nx,ny);
        }
        if (tilt) {
            iRSF q1("theta"), q2("phi");
            Q1.resize(nxyz); q1>>tt; expand3d(tt,Q1,nb,nz,nx,ny);
            Q2.resize(nxyz); q2>>tt; expand3d(tt,Q2,nb,nz,nx,ny);
        } else {
            Q1.resize(nxyz,0.);
            Q2.resize(nxyz,0.);
        }
    }
    /* convert degree to radian */
    for (int i=0; i<nxyz; i++) {
        Q1[i] *= SF_PI/180.;
        Q2[i] *= SF_PI/180.;
    }

    /*------------------------------------------------------------*/
    /* Fourier domain parameters                                  */
    /*------------------------------------------------------------*/
    /* Fourier spectra demension */
    int nkz,nkx,nky,nk;
    nkz=kiss_fft_next_fast_size(nzpad);
    nkx=kiss_fft_next_fast_size(nxpad);
    nky=kiss_fft_next_fast_size(nypad);
    nk = nky*nkx*nkz;

    float dkz,dkx,dky;
    float okz,okx,oky;
    dkz=2.*SF_PI/dz/nkz;
    dkx=2.*SF_PI/dx/nkx;
    dky=2.*SF_PI/dy/nky;
    okz=-SF_PI/dz;
    okx=-SF_PI/dx;
    oky=-SF_PI/dy;

    /* phase directions (unnormalized) */
    rkz.resize(nk);
    rkx.resize(nk);
    rky.resize(nk);
    for(int iy=0; iy<nky; iy++)
    {
        double ky = oky+iy*dky;
        if(ky==0.0) ky=SF_EPS*dky;

        for(int ix=0; ix<nkx; ix++)
        {
            double kx = okx+ix*dkx;
            if(kx==0.0) kx=SF_EPS*dkx;

            for(int iz=0; iz<nkz; iz++)
            {
                double kz = okz+iz*dkz;
                if(kz==0.0) kz=SF_EPS*dkz;

                int i = ((iy*nkx)+ix)*nkz+iz;
                rky[i] = ky;
                rkx[i] = kx;
                rkz[i] = kz;
            }
        }
    }

    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    if(verb) { 
        sf_warning("*******************************************");
        if(mode==0)      sf_warning("*** Elastic lowrank decomposition (P+S) ***");
        else if(mode==1) sf_warning("*** Elastic lowrank decomposition for P ***");
        else if(mode==2) sf_warning("*** Elastic lowrank decomposition for S ***");
        sf_warning("*******************************************");
        sf_warning("npk=%d,eps=%f",npk,eps);
        sf_warning("Reading velocity model parameters...");
        sf_warning("nz=%d, nx=%d, ny=%d",nz,nx,ny);
        sf_warning("dz=%f, dx=%f, dy=%f",dz,dx,dy);
        sf_warning("oz=%f, ox=%f, oy=%f",oz,ox,oy);
        sf_warning("nzpad=%d, nxpad=%d, nypad=%d",nzpad,nxpad,nypad);
        sf_warning("ozpad=%f, oxpad=%f, oypad=%f",ozpad,oxpad,oypad);
        sf_warning("nkz=%d, nkx=%d, nky=%d",nkz,nkx,nky);
        sf_warning("dkz=%f, dkx=%f, dky=%f",dkz,dkx,dky);
        sf_warning("okz=%f, okx=%f, oky=%f",okz,okx,oky);
        sf_warning("CPU time for preparation: %f(second)",timespent);
        sf_warning("*******************************************");
        sf_warning("*******************************************");
    }

    /*------------------------------------------------------------*/
    /* lowrank decomposition preparation                          */
    /*------------------------------------------------------------*/
    /* FIO dimension */
    int m=nxyz;   // number of rows (x)
    int n=nk;        // number of cols (k)
    int n2_sum=0;    // cumulative ranks
    int *n2s;
    n2s = sf_intalloc(6);

    /* col and row index arrays */
    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) midx[k] = k;
    for (int k=0; k < n; k++) nidx[k] = k;

    /* reference cols (left) and rows (right) index arrays */
    vector<int> lidx, ridx;
    ZpxNumMat mid;
    sf_complex *ldataxx, *ldatayy, *ldatazz, *ldataxy, *ldataxz, *ldatayz;
    sf_complex *rdataxx, *rdatayy, *rdatazz, *rdataxy, *rdataxz, *rdatayz;

    /* preparation for reduced lowrank */
    vector<int> ms, ns, js;
    ms.resize(3); ms[0] = nzpad; ms[1] = nxpad; ms[2] = nypad;
    ns.resize(3); ns[0] = nkz;   ns[1] = nkx;   ns[2] = nky;
    js.resize(3); js[0] = jump;  js[1] = jump;  js[2] = jump;

    /*------------------------------------------------------------*/
    /* lowrank decomposition and write to files                   */
    /*------------------------------------------------------------*/
    component = 0;
    /* xx component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,samplexx3,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of xx: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldataxx = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldataxx[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdataxx = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdataxx[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* yy component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of yy: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatayy = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatayy[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatayy = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatayy[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* zz component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of zz: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatazz = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatazz[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatazz = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatazz[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* xy component */
    {
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of xy: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldataxy = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldataxy[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdataxy = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdataxy[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* xz component */
    {
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of xz: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldataxz = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldataxz[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdataxz = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdataxz[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* yz component */
    {
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of yz: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatayz = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatayz[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatayz = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatayz[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }

    /*------------------------------------------------------------*/
    /* output results                                             */
    /*------------------------------------------------------------*/
    sf_file rank  = sf_output("out");
    sf_settype(rank,SF_INT);
    sf_putint(rank,"n1",6);
    sf_putint(rank,"n2",1);
    sf_putint(rank,"n3",1);
    sf_intwrite(n2s,6,rank);

    sf_file left = sf_output("left" );
    sf_settype(left,SF_COMPLEX);
    sf_putint(left,"n1",m);
    sf_putint(left,"n2",n2_sum);
    sf_putint(left,"n3",1);
    sf_complexwrite(ldataxx,m*n2s[0],left);
    sf_complexwrite(ldatayy,m*n2s[1],left);
    sf_complexwrite(ldatazz,m*n2s[2],left);
    sf_complexwrite(ldataxy,m*n2s[3],left);
    sf_complexwrite(ldataxz,m*n2s[4],left);
    sf_complexwrite(ldatayz,m*n2s[5],left);

    sf_file right = sf_output("right" );
    sf_settype(right,SF_COMPLEX);
    sf_putint(right,"n1",n);
    sf_putint(right,"n2",n2_sum);
    sf_putint(right,"n3",1);
    sf_complexwrite(rdataxx,n*n2s[0],right);
    sf_complexwrite(rdatayy,n*n2s[1],right);
    sf_complexwrite(rdatazz,n*n2s[2],right);
    sf_complexwrite(rdataxy,n*n2s[3],right);
    sf_complexwrite(rdataxz,n*n2s[4],right);
    sf_complexwrite(rdatayz,n*n2s[5],right);

    /*------------------------------------------------------------*/
    /* finalize calculation                                       */
    /*------------------------------------------------------------*/
    free(n2s);
    free(ldataxx); free(ldatayy); free(ldatazz); free(ldataxy); free(ldataxz); free(ldatayz);
    free(rdataxx); free(rdatayy); free(rdatazz); free(rdataxy); free(rdataxz); free(rdatayz);

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    if(verb) {
        sf_warning("CPU time for low-rank decomp: %f(second)",timespent);
        sf_warning("Total # of ranks: %d",n2_sum);
        sf_warning("lowrank decomposition completed, gracefully exiting...");
    }

    exit(0);
}

static int sample(vector<int>& rs, vector<int>& cs, ZpxNumMat& res)
/* component-wise 3d elastic wave extrapolation operator */
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);
    setvalue(res,zpx(0.0,0.0));

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double c11 = C11[i]; double c12 = C12[i]; double c13 = C13[i];
        double c22 = C22[i]; double c23 = C23[i]; double c33 = C33[i];
        double c44 = C44[i]; double c55 = C55[i]; double c66 = C66[i];
        double ss1 = sin(Q1[i]); double cc1 = cos(Q1[i]);
        double ss2 = sin(Q2[i]); double cc2 = cos(Q2[i]);
        double c14,c15,c16,c24,c25,c26,c34,c35,c36,c45,c46,c56;
        if (tric) {
            c14 = C14[i]; c15 = C15[i]; c16 = C16[i];
            c24 = C24[i]; c25 = C25[i]; c26 = C26[i];
            c34 = C34[i]; c35 = C35[i]; c36 = C36[i];
            c45 = C45[i]; c46 = C46[i]; c56 = C56[i];
        }

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            double kx  = kx0*cc2+ky0*ss2;
            double ky  =-kx0*ss2*cc1+ky0*cc2*cc1+kz0*ss1;
            double kz  = kx0*ss2*ss1-ky0*cc2*ss1+kz0*cc1;
            if(kx==0 && ky==0 && kz==0)
            {
                res(a,b) = zpx(0.0,0.0);
                continue;
            }

            double kx2 = kx*kx;
            double ky2 = ky*ky;
            double kz2 = kz*kz;
            double kxy = kx*ky;
            double kxz = kx*kz;
            double kyz = ky*kz;

            /* Christoffel matrix */
            if (tric) {
                Chr[0] =          c11*kx2 + c66*ky2 + c55*kz2 +     2*c16*kxy +     2*c15*kxz +     2*c56*kyz;
                Chr[4] =          c66*kx2 + c22*ky2 + c44*kz2 +     2*c26*kxy +     2*c46*kxz +     2*c24*kyz;
                Chr[8] =          c55*kx2 + c44*ky2 + c33*kz2 +     2*c45*kxy +     2*c35*kxz +     2*c34*kyz;
                Chr[1] = Chr[3] = c16*kx2 + c26*ky2 + c45*kz2 + (c12+c66)*kxy + (c14+c56)*kxz + (c25+c46)*kyz;
                Chr[2] = Chr[6] = c15*kx2 + c46*ky2 + c35*kz2 + (c14+c56)*kxy + (c13+c55)*kxz + (c36+c45)*kyz;
                Chr[5] = Chr[7] = c56*kx2 + c24*ky2 + c34*kz2 + (c25+c46)*kxy + (c36+c45)*kxz + (c23+c44)*kyz;
            } else {
                Chr[0] =          c11*kx2 + c66*ky2 + c55*kz2;
                Chr[4] =          c66*kx2 + c22*ky2 + c44*kz2;
                Chr[8] =          c55*kx2 + c44*ky2 + c33*kz2;
                Chr[1] = Chr[3] =                               (c12+c66)*kxy;
                Chr[2] = Chr[6] =                                               (c13+c55)*kxz;
                Chr[5] = Chr[7] =                                                               (c23+c44)*kyz;
            }

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            /* slow S wave */
            double v1t = sqrt(ww[0])*dt; // v_{s2}*k*dt
            double u1x = Chr[0];
            double u1y = Chr[1];
            double u1z = Chr[2];
            /* fast S wave */
            double v2t = sqrt(ww[1])*dt; // v_{s1}*k*dt
            double u2x = Chr[3];
            double u2y = Chr[4];
            double u2z = Chr[5];
            /* P wave */
            double v3t = sqrt(ww[2])*dt; // v_{p}*k*dt
            double u3x = Chr[6];
            double u3y = Chr[7];
            double u3z = Chr[8];
            //sf_warning("v1t=%f,u1x=%f,v2t=%f,u2x=%f,v3t=%f,u3x=%f",v1t,u1x,v2t,u2x,v3t,u3x);

            switch (component) {
                case 0: // xx
                    if(mode==0) {
                        if(tstp) {
                            if(pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3x,0) + zpx(-v2t*v2t*u2x*u2x,0) + zpx(-v1t*v1t*u1x*u1x,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3x,0) + zpx(2*(cos(v2t)-1.)*u2x*u2x,0) + zpx(2*(cos(v1t)-1.)*u1x*u1x,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3x + zpx(cos(v2t),sin(v2t))*u2x*u2x + zpx(cos(v1t),sin(v1t))*u1x*u1x ;
                    } else if(mode==1) {
                        if(tstp) {
                            if(pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3x,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3x,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3x ;
                    } else if(mode==2) {
                        if(tstp) {
                            if(pseu) res(a,b) = zpx(-v2t*v2t*u2x*u2x,0) + zpx(-v1t*v1t*u1x*u1x,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2x*u2x,0) + zpx(2*(cos(v1t)-1.)*u1x*u1x,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2x*u2x + zpx(cos(v1t),sin(v1t))*u1x*u1x ;
                    }
                    break;
                case 1: // yy
                    if(mode==0) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3y*u3y,0) + zpx(-v2t*v2t*u2y*u2y,0) + zpx(-v1t*v1t*u1y*u1y,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3y*u3y,0) + zpx(2*(cos(v2t)-1.)*u2y*u2y,0) + zpx(2*(cos(v1t)-1.)*u1y*u1y,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3y*u3y + zpx(cos(v2t),sin(v2t))*u2y*u2y + zpx(cos(v1t),sin(v1t))*u1y*u1y ;
                    } else if(mode==1) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3y*u3y,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3y*u3y,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3y*u3y ;
                    } else if(mode==2) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v2t*v2t*u2y*u2y,0) + zpx(-v1t*v1t*u1y*u1y,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2y*u2y,0) + zpx(2*(cos(v1t)-1.)*u1y*u1y,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2y*u2y + zpx(cos(v1t),sin(v1t))*u1y*u1y ;
                    }
                    break;
                case 2: // zz
                    if(mode==0) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3z*u3z,0) + zpx(-v2t*v2t*u2z*u2z,0) + zpx(-v1t*v1t*u1z*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3z*u3z,0) + zpx(2*(cos(v2t)-1.)*u2z*u2z,0) + zpx(2*(cos(v1t)-1.)*u1z*u1z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3z*u3z + zpx(cos(v2t),sin(v2t))*u2z*u2z + zpx(cos(v1t),sin(v1t))*u1z*u1z ;
                    } else if(mode==1) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3z*u3z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3z*u3z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3z*u3z ;
                    } else if(mode==2) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v2t*v2t*u2z*u2z,0) + zpx(-v1t*v1t*u1z*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2z*u2z,0) + zpx(2*(cos(v1t)-1.)*u1z*u1z,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2z*u2z + zpx(cos(v1t),sin(v1t))*u1z*u1z ;
                    }
                    break;
                case 3: //xy
                    if(mode==0) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3y,0) + zpx(-v2t*v2t*u2x*u2y,0) + zpx(-v1t*v1t*u1x*u1y,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3y,0) + zpx(2*(cos(v2t)-1.)*u2x*u2y,0) + zpx(2*(cos(v1t)-1.)*u1x*u1y,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3y + zpx(cos(v2t),sin(v2t))*u2x*u2y + zpx(cos(v1t),sin(v1t))*u1x*u1y ;
                    } else if(mode==1) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3y,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3y,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3y ;
                    } else if(mode==2) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v2t*v2t*u2x*u2y,0) + zpx(-v1t*v1t*u1x*u1y,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2x*u2y,0) + zpx(2*(cos(v1t)-1.)*u1x*u1y,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2x*u2y + zpx(cos(v1t),sin(v1t))*u1x*u1y ;
                    }
                    break;
                case 4: // xz
                    if(mode==0) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3z,0) + zpx(-v2t*v2t*u2x*u2z,0) + zpx(-v1t*v1t*u1x*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3z,0) + zpx(2*(cos(v2t)-1.)*u2x*u2z,0) + zpx(2*(cos(v1t)-1.)*u1x*u1z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3z + zpx(cos(v2t),sin(v2t))*u2x*u2z + zpx(cos(v1t),sin(v1t))*u1x*u1z ;
                    } else if(mode==1) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3x*u3z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3x*u3z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3x*u3z ;
                    } else if(mode==2) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v2t*v2t*u2x*u2z,0) + zpx(-v1t*v1t*u1x*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2x*u2z,0) + zpx(2*(cos(v1t)-1.)*u1x*u1z,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2x*u2z + zpx(cos(v1t),sin(v1t))*u1x*u1z ;
                    }
                    break;
                case 5: // yz
                    if(mode==0) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3y*u3z,0) + zpx(-v2t*v2t*u2y*u2z,0) + zpx(-v1t*v1t*u1y*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3y*u3z,0) + zpx(2*(cos(v2t)-1.)*u2y*u2z,0) + zpx(2*(cos(v1t)-1.)*u1y*u1z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3y*u3z + zpx(cos(v2t),sin(v2t))*u2y*u2z + zpx(cos(v1t),sin(v1t))*u1y*u1z ;
                    } else if(mode==1) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v3t*v3t*u3y*u3z,0) ;
                            else res(a,b) = zpx(2*(cos(v3t)-1.)*u3y*u3z,0) ;
                        } else res(a,b) = zpx(cos(v3t),sin(v3t))*u3y*u3z ;
                    } else if(mode==2) {
                        if (tstp) {
                            if (pseu) res(a,b) = zpx(-v2t*v2t*u2y*u2z,0) + zpx(-v1t*v1t*u1y*u1z,0) ;
                            else res(a,b) = zpx(2*(cos(v2t)-1.)*u2y*u2z,0) + zpx(2*(cos(v1t)-1.)*u1y*u1z,0) ;
                        } else res(a,b) = zpx(cos(v2t),sin(v2t))*u2y*u2z + zpx(cos(v1t),sin(v1t))*u1y*u1z ;
                    }
                    break;
                default:
                    break;
            }

        }// b loop
    }// a loop

    return 0;
}

void expand3d(std::valarray<float>  a, 
	      std::valarray<float>& b, 
              int nb,
              int nz,
	      int nx,
              int ny)
/*< expand domain >*/
{
    int iz,ix,iy;
    int nzpad = nz+2*nb;
    int nxpad = nx+2*nb;
    int nypad = ny+2*nb;

    for         (iy=0;iy<ny;iy++) {
	for     (ix=0;ix<nx;ix++) {
	    for (iz=0;iz<nz;iz++) {
		//b[nb+iy][nb+ix][nb+iz] = a[iy][ix][iz];
                b[((nb+iy)*nxpad+nb+ix)*nzpad+nb+iz] = a[(iy*nx+ix)*nz+iz];
	    }
	}
    }

    for         (iy=0; iy<nypad; iy++) {
	for     (ix=0; ix<nxpad; ix++) {
	    for (iz=0; iz<nb;    iz++) {
		//b[iy][ix][      iz  ] = b[iy][ix][      nb  ];
		//b[iy][ix][nzpad-iz-1] = b[iy][ix][nzpad-nb-1];
                b[(iy*nxpad+ix)*nzpad+      iz  ] = b[(iy*nxpad+ix)*nzpad+      nb  ];
                b[(iy*nxpad+ix)*nzpad+nzpad-iz-1] = b[(iy*nxpad+ix)*nzpad+nzpad-nb-1];
	    }
	}
    }


    for         (iy=0; iy<nypad; iy++) {
	for     (ix=0; ix<nb;    ix++) {
	    for (iz=0; iz<nzpad; iz++) {
		//b[iy][      ix  ][iz] = b[iy][      nb  ][iz];
		//b[iy][nxpad-ix-1][iz] = b[iy][nxpad-nb-1][iz];
                b[(iy*nxpad+      ix  )*nzpad+iz] = b[(iy*nxpad+      nb  )*nzpad+iz];
                b[(iy*nxpad+nxpad-ix-1)*nzpad+iz] = b[(iy*nxpad+nxpad-nb-1)*nzpad+iz];
	    }
	}
    }

    for         (iy=0; iy<nb;    iy++) {
	for     (ix=0; ix<nxpad; ix++) {
	    for (iz=0; iz<nzpad; iz++) {
		//b[      iy  ][ix][iz] = b[      nb  ][ix][iz];
		//b[nypad-iy-1][ix][iz] = b[nypad-nb-1][ix][iz];
                b[((      iy  )*nxpad+ix)*nzpad+iz] = b[((      nb  )*nxpad+ix)*nzpad+iz];
                b[((nypad-iy-1)*nxpad+ix)*nzpad+iz] = b[((nypad-nb-1)*nxpad+ix)*nzpad+iz];
	    }
	}
    }

}

