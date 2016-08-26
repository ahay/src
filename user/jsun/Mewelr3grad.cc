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

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;
// not a namespace name?
//using namespace std::complex_literals;

/* automatically generated header files */
extern "C"{
#include "eigen3x3.h"
}

/* centered derivatives */
#define NOP 4 /* derivative operator half-size */
#define C1 4.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 4.0f/105.0f
#define C4 -1.0f/280.0f
#define Dy(a,iy,ix,iz,nx,nz,s) (C4*(a[((iy+4)*nx+ix)*nz+iz] - a[((iy-4)*nx+ix)*nz+iz]) +    \
                                C3*(a[((iy+3)*nx+ix)*nz+iz] - a[((iy-3)*nx+ix)*nz+iz]) +    \
                                C2*(a[((iy+2)*nx+ix)*nz+iz] - a[((iy-2)*nx+ix)*nz+iz]) +    \
                                C1*(a[((iy+1)*nx+ix)*nz+iz] - a[((iy-1)*nx+ix)*nz+iz])  )*s
#define Dx(a,iy,ix,iz,nx,nz,s) (C4*(a[(iy*nx+(ix+4))*nz+iz] - a[(iy*nx+(ix-4))*nz+iz]) +    \
                                C3*(a[(iy*nx+(ix+3))*nz+iz] - a[(iy*nx+(ix-3))*nz+iz]) +    \
                                C2*(a[(iy*nx+(ix+2))*nz+iz] - a[(iy*nx+(ix-2))*nz+iz]) +    \
                                C1*(a[(iy*nx+(ix+1))*nz+iz] - a[(iy*nx+(ix-1))*nz+iz])  )*s
#define Dz(a,iy,ix,iz,nx,nz,s) (C4*(a[(iy*nx+ix)*nz+(iz+4)] - a[(iy*nx+ix)*nz+(iz-4)]) +    \
                                C3*(a[(iy*nx+ix)*nz+(iz+3)] - a[(iy*nx+ix)*nz+(iz-3)]) +    \
                                C2*(a[(iy*nx+ix)*nz+(iz+2)] - a[(iy*nx+ix)*nz+(iz-2)]) +    \
                                C1*(a[(iy*nx+ix)*nz+(iz+1)] - a[(iy*nx+ix)*nz+(iz-1)])  )*s

/* static variables */
static float dt;
static bool tric,tstp,pseu,grad;
static std::valarray<float>  C11,C12,C13,C22,C23,C33,C44,C55,C66,Q1,Q2;
static std::valarray<float>  C11dx,C12dx,C12dy,C13dx,C13dz,C22dy,C23dy,C23dz,C33dz,C44dy,C44dz,C55dx,C55dz,C66dx,C66dy;
static std::valarray<float>  C14,C15,C16,C24,C25,C26,C34,C35,C36,C45,C46,C56;
static std::valarray<double> rkz, rkx, rky;
static int component, mode;

/* subroutines */
static int sample(vector<int>& rs, vector<int>& cs, ZpxNumMat& res);

void expand3d(std::valarray<float>  a, std::valarray<float>& b, int nb, int nz, int nx, int ny);
void dr3d(std::valarray<float>  u, std::valarray<float>& du, int nz, int nx, int ny, float dz, float dx, float dy, int which);


/* global variables needed by lapack SVD routines */
char    jobz='V';
char    uplo='U';
//int     M=3;
//int     LDA=M;
//int     LWORK=8*M;
//int     INFO;
double  Chr[9], ww[3], work[24];

char    JOBVL='N';
char    JOBVR='V';
int     M=3;
int     LDA=M;
int     LDVL=3;
int     LDVR=3;
int     LWORK=8*M;
int     INFO;
zpx     Chr2[9], ww2[3], work2[24], vl[9], vr[9];
double  rwork[6];

int     ipiv[3];

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
    par.get("grad",grad,false);      // include gradient term
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
    /* calculate Cij gradient */
    C11dx.resize(nxyz); dr3d(C11, C11dx, nzpad, nxpad, nypad, dz, dx, dy, 2);
    C12dx.resize(nxyz); dr3d(C12, C12dx, nzpad, nxpad, nypad, dz, dx, dy, 2);
    C12dy.resize(nxyz); dr3d(C12, C12dy, nzpad, nxpad, nypad, dz, dx, dy, 3);
    C13dx.resize(nxyz); dr3d(C13, C13dx, nzpad, nxpad, nypad, dz, dx, dy, 2);
    C13dz.resize(nxyz); dr3d(C13, C13dz, nzpad, nxpad, nypad, dz, dx, dy, 1);
    C22dy.resize(nxyz); dr3d(C22, C22dy, nzpad, nxpad, nypad, dz, dx, dy, 3);
    C23dy.resize(nxyz); dr3d(C23, C23dy, nzpad, nxpad, nypad, dz, dx, dy, 3);
    C23dz.resize(nxyz); dr3d(C23, C23dz, nzpad, nxpad, nypad, dz, dx, dy, 1);
    C33dz.resize(nxyz); dr3d(C33, C33dz, nzpad, nxpad, nypad, dz, dx, dy, 1);
    C44dy.resize(nxyz); dr3d(C44, C44dy, nzpad, nxpad, nypad, dz, dx, dy, 3);
    C44dz.resize(nxyz); dr3d(C44, C44dz, nzpad, nxpad, nypad, dz, dx, dy, 1);
    C55dx.resize(nxyz); dr3d(C55, C55dx, nzpad, nxpad, nypad, dz, dx, dy, 2);
    C55dz.resize(nxyz); dr3d(C55, C55dz, nzpad, nxpad, nypad, dz, dx, dy, 1);
    C66dx.resize(nxyz); dr3d(C66, C66dx, nzpad, nxpad, nypad, dz, dx, dy, 2);
    C66dy.resize(nxyz); dr3d(C66, C66dy, nzpad, nxpad, nypad, dz, dx, dy, 3);
    //oRSF fC11dx("C11dx"); fC11dx.put("n1",nzpad); fC11dx.put("n2",nxpad); fC11dx.put("n3",nypad); fC11dx<<C11dx;
    //oRSF fC12dy("C12dy"); fC12dy.put("n1",nzpad); fC12dy.put("n2",nxpad); fC12dy.put("n3",nypad); fC12dy<<C12dy;
    //oRSF fC13dz("C13dz"); fC13dz.put("n1",nzpad); fC13dz.put("n2",nxpad); fC13dz.put("n3",nypad); fC13dz<<C13dz;
    //exit(0);
 
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
        sf_warning("****** Elastic lowrank decomposition ******");
        sf_warning("*******************************************");
        if (grad) sf_warning("Lowrank for RITE including gradient term!");
        else sf_warning("Lowrank for RITE without gradient term!");
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
    n2s = sf_intalloc(9);

    /* col and row index arrays */
    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) midx[k] = k;
    for (int k=0; k < n; k++) nidx[k] = k;

    /* reference cols (left) and rows (right) index arrays */
    vector<int> lidx, ridx;
    ZpxNumMat mid;
    sf_complex *ldataxx, *ldataxy, *ldataxz, *ldatayx, *ldatayy, *ldatayz, *ldatazx, *ldatazy, *ldatazz;
    sf_complex *rdataxx, *rdataxy, *rdataxz, *rdatayx, *rdatayy, *rdatayz, *rdatazx, *rdatazy, *rdatazz;

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
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
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
    /* xy component */
    {
        srand48(seed);
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
        srand48(seed);
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
    /* yx component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of yx: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatayx = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatayx[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatayx = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatayx[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
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
    /* yz component */
    {
        srand48(seed);
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
    /* zx component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of zx: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatazx = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatazx[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatazx = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatazx[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
            }
        }

        n2s[component] = n2;
        n2_sum += n2;
        component++;
    }
    /* zy component */
    {
        srand48(seed);
        //iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );
        iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );

        int m2=mid.m();
        int n2=mid.n();
        if(verb) sf_warning("Rank of zy: # of ref. cols(k)=%d; # of ref. rows(x)=%d",m2,n2);

        /* left */
        ZpxNumMat lmat(m,m2);
        iC ( sample(midx,lidx,lmat) );

        ZpxNumMat lmat2(m,n2);
        iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

        zpx *ldat = lmat2.data();
        ldatazy = sf_complexalloc(m*n2);
        for (int k=0; k < m*n2; k++) 
            ldatazy[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

        /* right */
        ZpxNumMat rmat(n2,n);
        iC ( sample(ridx,nidx,rmat) );

        zpx *rdat = rmat.data();
        rdatazy = sf_complexalloc(n*n2);
        for     (int i=0; i < n;  i++) {
            for (int j=0; j < n2; j++) {
                rdatazy[j*n+i] = sf_cmplx(real(rdat[i*n2+j]),imag(rdat[i*n2+j])); // transposition
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

    /*------------------------------------------------------------*/
    /* output results                                             */
    /*------------------------------------------------------------*/
    sf_file rank  = sf_output("out");
    sf_settype(rank,SF_INT);
    sf_putint(rank,"n1",9);
    sf_putint(rank,"n2",1);
    sf_putint(rank,"n3",1);
    sf_intwrite(n2s,9,rank);

    sf_file left = sf_output("left" );
    sf_settype(left,SF_COMPLEX);
    sf_putint(left,"n1",m);
    sf_putint(left,"n2",n2_sum);
    sf_putint(left,"n3",1);
    sf_complexwrite(ldataxx,m*n2s[0],left);
    sf_complexwrite(ldataxy,m*n2s[1],left);
    sf_complexwrite(ldataxz,m*n2s[2],left);
    sf_complexwrite(ldatayx,m*n2s[3],left);
    sf_complexwrite(ldatayy,m*n2s[4],left);
    sf_complexwrite(ldatayz,m*n2s[5],left);
    sf_complexwrite(ldatazx,m*n2s[6],left);
    sf_complexwrite(ldatazy,m*n2s[7],left);
    sf_complexwrite(ldatazz,m*n2s[8],left);

    sf_file right = sf_output("right" );
    sf_settype(right,SF_COMPLEX);
    sf_putint(right,"n1",n);
    sf_putint(right,"n2",n2_sum);
    sf_putint(right,"n3",1);
    sf_complexwrite(rdataxx,n*n2s[0],right);
    sf_complexwrite(rdataxy,n*n2s[1],right);
    sf_complexwrite(rdataxz,n*n2s[2],right);
    sf_complexwrite(rdatayx,n*n2s[3],right);
    sf_complexwrite(rdatayy,n*n2s[4],right);
    sf_complexwrite(rdatayz,n*n2s[5],right);
    sf_complexwrite(rdatazx,n*n2s[6],right);
    sf_complexwrite(rdatazy,n*n2s[7],right);
    sf_complexwrite(rdatazz,n*n2s[8],right);

    /*------------------------------------------------------------*/
    /* finalize calculation                                       */
    /*------------------------------------------------------------*/
    free(n2s);
    free(ldataxx); free(ldataxy); free(ldataxz); free(ldatayx); free(ldatayy); free(ldatayz); free(ldatazx); free(ldatazy); free(ldatazz);
    free(rdataxx); free(rdataxy); free(rdataxz); free(rdatayx); free(rdatayy); free(rdatayz); free(rdatazx); free(rdatazy); free(rdatazz);

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
        double c11dx = C11dx[i]; double c12dx = C12dx[i]; double c12dy = C12dy[i]; 
        double c13dx = C13dx[i]; double c13dz = C13dz[i]; double c22dy = C22dy[i]; 
        double c23dy = C23dy[i]; double c23dz = C23dz[i]; double c33dz = C33dz[i]; 
        double c44dy = C44dy[i]; double c44dz = C44dz[i]; double c55dx = C55dx[i]; 
        double c55dz = C55dz[i]; double c66dx = C66dx[i]; double c66dy = C66dy[i];

        //c11dx=c12dx=c12dy=c13dx=c13dz=c22dy=c23dy=c23dz=c33dz=c44dy=c44dz=c55dx=c55dz=c66dx=c66dy=20;

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

            zpx v1t,u1x,u1y,u1z,iu1x,iu1y,iu1z;
            zpx v2t,u2x,u2y,u2z,iu2x,iu2y,iu2z;
            zpx v3t,u3x,u3y,u3z,iu3x,iu3y,iu3z;

            if (grad) {
                Chr2[0] = zpx(Chr[0],-(c11dx*kx + c66dy*ky + c55dz*kz)); Chr2[1] = zpx(Chr[1],-(c12dx*ky + c66dy*kx)); Chr2[2] = zpx(Chr[2],-(c13dx*kz + c55dz*kx));
                Chr2[3] = zpx(Chr[3],-(c12dy*kx + c66dx*ky)); Chr2[4] = zpx(Chr[4],-(c66dx*kx + c22dy*ky + c44dz*kz)); Chr2[5] = zpx(Chr[5],-(c23dy*kz + c44dz*ky));
                Chr2[6] = zpx(Chr[6],-(c13dz*kx + c55dx*kz)); Chr2[7] = zpx(Chr[7],-(c23dz*ky + c44dy*kz)); Chr2[8] = zpx(Chr[8],-(c55dx*kx + c44dy*ky + c33dz*kz));
                /*
                Chr2[0] = zpx(Chr[0],0); Chr2[1] = zpx(Chr[1],0); Chr2[2] = zpx(Chr[2],0);
                Chr2[3] = zpx(Chr[3],0); Chr2[4] = zpx(Chr[4],0); Chr2[5] = zpx(Chr[5],0);
                Chr2[6] = zpx(Chr[6],0); Chr2[7] = zpx(Chr[7],0); Chr2[8] = zpx(Chr[8],0);
                */

                // LAPACK's eigen-decom (complex version)
                zgeev_(&JOBVL, &JOBVR, &M, (MKL_Complex16*) Chr2, &LDA, (MKL_Complex16*) ww2,  (MKL_Complex16*) vl, &LDVL,  (MKL_Complex16*) vr, &LDVR,  (MKL_Complex16*) work2, &LWORK, rwork, &INFO);

                zpx ivr[9];
                for (int id=0; id<9; id++) ivr[id]=vr[id];
                zgetrf_(&M, &M,  (MKL_Complex16*) ivr, &LDA, ipiv, &INFO);
                zgetri_(&M,  (MKL_Complex16*) ivr, &LDA, ipiv,  (MKL_Complex16*) work2, &LWORK, &INFO);

                /* slow S wave */
                v1t = sqrt(ww2[0])*zpx(dt,0); // v_{s2}*k*dt
                u1x = vr[0]; iu1x = ivr[0];
                u1y = vr[1]; iu1y = ivr[3];
                u1z = vr[2]; iu1z = ivr[6];
                /* fast S wave */
                v2t = sqrt(ww2[1])*zpx(dt,0); // v_{s1}*k*dt
                u2x = vr[3]; iu2x = ivr[1];
                u2y = vr[4]; iu2y = ivr[4];
                u2z = vr[5]; iu2z = ivr[7];
                /* P wave */
                v3t = sqrt(ww2[2])*zpx(dt,0); // v_{p}*k*dt
                u3x = vr[6]; iu3x = ivr[2];
                u3y = vr[7]; iu3y = ivr[5];
                u3z = vr[8]; iu3z = ivr[8];

                if (real(v1t)<0 || real(v2t)<0 || real(v3t)<0) sf_error("Negative velocity detected: real(v1t)=%f, real(v2t)=%f, real(v3t)=%f.",real(v1t),real(v2t),real(v3t));

            } else {
                // LAPACK's eigen-decom (real-version)
                dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
                /* slow S wave */
                v1t = sqrt(ww[0])*zpx(dt,0); // v_{s2}*k*dt
                u1x = Chr[0]; iu1x = Chr[0];
                u1y = Chr[1]; iu1y = Chr[1];
                u1z = Chr[2]; iu1z = Chr[2];
                /* fast S wave */
                v2t = sqrt(ww[1])*zpx(dt,0); // v_{s1}*k*dt
                u2x = Chr[3]; iu2x = Chr[3];
                u2y = Chr[4]; iu2y = Chr[4];
                u2z = Chr[5]; iu2z = Chr[5];
                /* P wave */
                v3t = sqrt(ww[2])*zpx(dt,0); // v_{p}*k*dt
                u3x = Chr[6]; iu3x = Chr[6];
                u3y = Chr[7]; iu3y = Chr[7];
                u3z = Chr[8]; iu3z = Chr[8];
            }
            
            // test inverse
            //sf_warning("Row1: (%f,%f) (%f,%f) (%f,%f)",real(u1x*iu1x+u1y*iu1y+u1z*iu1z),imag(u1x*iu1x+u1y*iu1y+u1z*iu1z),real(u1x*iu2x+u1y*iu2y+u1z*iu2z),imag(u1x*iu2x+u1y*iu2y+u1z*iu2z),real(u1x*iu3x+u1y*iu3y+u1z*iu3z),imag(u1x*iu3x+u1y*iu3y+u1z*iu3z));
            //sf_warning("Row2: (%f,%f) (%f,%f) (%f,%f)",real(u2x*iu1x+u2y*iu1y+u2z*iu1z),imag(u2x*iu1x+u2y*iu1y+u2z*iu1z),real(u2x*iu2x+u2y*iu2y+u2z*iu2z),imag(u2x*iu2x+u2y*iu2y+u2z*iu2z),real(u2x*iu3x+u2y*iu3y+u2z*iu3z),imag(u2x*iu3x+u2y*iu3y+u2z*iu3z));
            //sf_warning("Row3: (%f,%f) (%f,%f) (%f,%f)",real(u3x*iu1x+u3y*iu1y+u3z*iu1z),imag(u3x*iu1x+u3y*iu1y+u3z*iu1z),real(u3x*iu2x+u3y*iu2y+u3z*iu2z),imag(u3x*iu2x+u3y*iu2y+u3z*iu2z),real(u3x*iu3x+u3y*iu3y+u3z*iu3z),imag(u3x*iu3x+u3y*iu3y+u3z*iu3z));
            // display eigenvectors
            //sf_warning("for 1: (%f,%f) (%f,%f) (%f,%f)",real(u1x),imag(u1x),real(u1y),imag(u1y),real(u1z),imag(u1z));
            //sf_warning("inv 1: (%f,%f) (%f,%f) (%f,%f)",real(iu1x),imag(iu1x),real(iu1y),imag(iu1y),real(iu1z),imag(iu1z));
            //sf_warning("for 2: (%f,%f) (%f,%f) (%f,%f)",real(u2x),imag(u2x),real(u2y),imag(u2y),real(u2z),imag(u2z));
            //sf_warning("inv 2: (%f,%f) (%f,%f) (%f,%f)",real(iu2x),imag(iu2x),real(iu2y),imag(iu2y),real(iu2z),imag(iu2z));
            //sf_warning("for 3: (%f,%f) (%f,%f) (%f,%f)",real(u3x),imag(u3x),real(u3y),imag(u3y),real(u3z),imag(u3z));
            //sf_warning("inv 3: (%f,%f) (%f,%f) (%f,%f)",real(iu3x),imag(iu3x),real(iu3y),imag(iu3y),real(iu3z),imag(iu3z));
            // sanity check using old eigen decomposition
            //sf_warning("1:v1t=%f,u1x=%f,v2t=%f,u2x=%f,v3t=%f,u3x=%f",sqrt(ww[0])*dt,Chr[0],sqrt(ww[1])*dt,Chr[3],sqrt(ww[2])*dt,Chr[6]);
            //sf_warning("2:v1t=%f,u1x=%f,v2t=%f,u2x=%f,v3t=%f,u3x=%f",real(sqrt(ww2[0])*zpx(dt,0)),real(u1x),real(sqrt(ww2[1])*zpx(dt,0)),real(u2x),real(sqrt(ww2[2])*zpx(dt,0)),real(u3x));
            // why the following won't work?
            //std::cout << "v1t = " << v1t << '\n';
            //std::cout << "v2t = " << v2t << '\n';
            //std::cout << "v3t = " << v3t << '\n';

            switch (component) {
                case 0: // xx
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3x + exp(v2t*zpx(0.,1.))*u2x*iu2x + exp(v1t*zpx(0.,1.))*u1x*iu1x;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3x                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2x*iu2x + exp(v1t*zpx(0.,1.))*u1x*iu1x;
                    break;
                case 1: // xy
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3y + exp(v2t*zpx(0.,1.))*u2x*iu2y + exp(v1t*zpx(0.,1.))*u1x*iu1y;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3y                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2x*iu2y + exp(v1t*zpx(0.,1.))*u1x*iu1y;
                    break;
                case 2: // xz
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3z + exp(v2t*zpx(0.,1.))*u2x*iu2z + exp(v1t*zpx(0.,1.))*u1x*iu1z;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3x*iu3z                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2x*iu2z + exp(v1t*zpx(0.,1.))*u1x*iu1z;
                    break;
                case 3: // yx
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3x + exp(v2t*zpx(0.,1.))*u2y*iu2x + exp(v1t*zpx(0.,1.))*u1y*iu1x;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3x                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2y*iu2x + exp(v1t*zpx(0.,1.))*u1y*iu1x;
                    break;
                case 4: // yy
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3y + exp(v2t*zpx(0.,1.))*u2y*iu2y + exp(v1t*zpx(0.,1.))*u1y*iu1y;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3y                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2y*iu2y + exp(v1t*zpx(0.,1.))*u1y*iu1y;
                    break;
                case 5: // yz
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3z + exp(v2t*zpx(0.,1.))*u2y*iu2z + exp(v1t*zpx(0.,1.))*u1y*iu1z;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3y*iu3z                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2y*iu2z + exp(v1t*zpx(0.,1.))*u1y*iu1z;
                    break;
                case 6: // zx
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3x + exp(v2t*zpx(0.,1.))*u2z*iu2x + exp(v1t*zpx(0.,1.))*u1z*iu1x;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3x                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2z*iu2x + exp(v1t*zpx(0.,1.))*u1z*iu1x;
                    break;
                case 7: // zy
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3y + exp(v2t*zpx(0.,1.))*u2z*iu2y + exp(v1t*zpx(0.,1.))*u1z*iu1y;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3y                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2z*iu2y + exp(v1t*zpx(0.,1.))*u1z*iu1y;
                    break;
                case 8: // zz
                    if     (mode==0) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3z + exp(v2t*zpx(0.,1.))*u2z*iu2z + exp(v1t*zpx(0.,1.))*u1z*iu1z;
                    else if(mode==1) res(a,b) = exp(v3t*zpx(0.,1.))*u3z*iu3z                                                              ;
                    else if(mode==2) res(a,b) =                                exp(v2t*zpx(0.,1.))*u2z*iu2z + exp(v1t*zpx(0.,1.))*u1z*iu1z;
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
                b[((nb+iy)*nxpad+nb+ix)*nzpad+nb+iz] = a[(iy*nx+ix)*nz+nz];
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

void dr3d(std::valarray<float>  u, 
          std::valarray<float>& du,
          int nz,
          int nx,
          int ny,
          float dz,
          float dx,
          float dy,
          int which)
/*< First order derivative in Z >*/
{
    int iz,ix,iy;
    float idd;

    switch (which) {
        case 2:
            idd = 1./dx;
            for         (iy=NOP;iy<ny-NOP;iy++) {
                for     (ix=NOP;ix<nx-NOP;ix++) {
                    for (iz=NOP;iz<nz-NOP;iz++) {
                        du[(iy*nx+ix)*nz+iz] = Dx(u,iy,ix,iz,nx,nz,idd);
                    }
                }
            }
            break;
        case 3:
            idd = 1./dy;
            for         (iy=NOP;iy<ny-NOP;iy++) {
                for     (ix=NOP;ix<nx-NOP;ix++) {
                    for (iz=NOP;iz<nz-NOP;iz++) {
                        du[(iy*nx+ix)*nz+iz] = Dy(u,iy,ix,iz,nx,nz,idd);
                    }
                }
            }
            break;
        default:
            idd = 1./dz;
            for         (iy=NOP;iy<ny-NOP;iy++) {
                for     (ix=NOP;ix<nx-NOP;ix++) {
                    for (iz=NOP;iz<nz-NOP;iz++) {
                        du[(iy*nx+ix)*nz+iz] = Dz(u,iy,ix,iz,nx,nz,idd);
                    }
                }
            }
    }

}
