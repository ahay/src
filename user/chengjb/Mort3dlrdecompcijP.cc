/* 3-D three-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and vector decomposition using low-rank symbol approximation in 3D ORT media (qP).

   Authors: Jiubing Cheng (Tongji University) and Sergey Fomel (The University of Texas at Austin)
   Modified: Yanadet Sripanich (The University of Texas at Austin)
     
   Copyright (C) 2012 Tongji University, Shanghai, China 
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
#include <rsf.hh>
#include <rsf.h>
#include <assert.h>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

/* prepared head files by myself */
#include "_cjb.h"

/* head files automatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h" // CHANGE IF MOVE TO MY FOLDER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "kykxkztaper.h"
#include "eigen3x3.h"
#include "puthead.h"
#include "decomplowrank.h"
}

static std::valarray<float> c11, c22, c33, c12, c13, c23, c44, c55, c66, phx, phy, phz;

static std::valarray<double> rkx, rky, rkz, rkk;

/* for low-rank decomp. */
static int samplebpxx(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebpyy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebpzz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebpxy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebpxz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebpyz(vector<int>& rs, vector<int>& cs, DblNumMat& res);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

/* definition for LAPACK SVD ROUTINEs */
char    jobz='V';  // for SVD 
char    uplo='U';  // for SVD 
int     M=3;       // for SVD 
int     LDA=M;     // for SVD 
int     LWORK=8*M; // for SVD 
int     INFO;      // for SVD 
double  Chr[9], ww[9], work[24];  // Lapack SVD array 

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4, t5, t44;
   float   timespent;

   t1=clock();

   iRSF par(0);
   int seed;
   par.get("seed",seed,time(NULL)); // seed for random number generator
   srand48(seed);

   float eps;
   par.get("eps",eps,1.e-6); // tolerance
       
   int npk;
   par.get("npk",npk,20); // maximum rank

   int ireconstruct;   // flag for reconstruct the W matrix or not

   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF c11r, c22r("c22"), c33r("c33"), c12r("c12"), c13r("c13"),c23r("c23"), c44r("c44"), c55r("c55"), c66r("c66"), phix("phix"), phiz("phiz"), phiy("phiy");

   /* Read/Write axes */
   int nxv, nyv, nzv;
   c11r.get("n1",nzv);
   c11r.get("n2",nxv);
   c11r.get("n3",nyv);

   float a1, a2, a3;
   c11r.get("o1",a1);
   c11r.get("o2",a2);
   c11r.get("o3",a3);

   float fxv, fyv, fzv;
   fxv=a2*1000.0;
   fyv=a3*1000.0;
   fzv=a1*1000.0;

   float dxv, dyv, dzv;
   c11r.get("d1",a1);
   c11r.get("d2",a2);
   c11r.get("d3",a3);
   dzv = a1*1000.0;
   dxv = a2*1000.0;
   dyv = a3*1000.0;

   /* Read/Write axes from Wavefields*/
   sf_file Fx, Fy, Fz;
   sf_axis az, ax, ay;
   
   Fx = sf_input("Elasticx");
   Fy = sf_input("Elasticy");
   Fz = sf_input("Elasticz");

   int   nx, ny, nz;
   float fx, fy, fz;
   float dx, dy, dz;

   az = sf_iaxa(Fx,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
   ax = sf_iaxa(Fx,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
   ay = sf_iaxa(Fx,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
   fy = sf_o(ay)*1000.0;
   fx = sf_o(ax)*1000.0;
   fz = sf_o(az)*1000.0;

   if(nx!=nxv || ny!=nyv || nz!=nzv){
     sf_warning("Dimension not match between model and data !");
     sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
     sf_warning("nxv=%d nyv=%d nzv=%d ",nxv,nyv,nzv);
     exit(0);
   }
   if(fabs(fx-fxv)>0.1 || fabs(fy-fyv)>0.1 || fabs(fz-fzv)>0.1){
     sf_warning("Coorinate original point not match between model and data !");
     sf_warning("fx=%d fy=%d fz=%d ",fx,fy,fz);
     sf_warning("fxv=%d fyv=%d fzv=%d ",fxv,fyv,fzv);
   }
   if(fabs(dx-dxv)>0.1 || fabs(dy-dyv)>0.1 || fabs(dz-dzv)>0.1){
     sf_warning("Sampling step not match between model and data !");
     sf_warning("dx=%d dy=%d dz=%d ",dx,dy,dz);
     sf_warning("dxv=%d dyv=%d dzv=%d ",dxv,dyv,dzv);
   }

   sf_warning("fx=%f fy=%f fz=%f dx=%f dy=%f dz=%f",fx,fy,fz,dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d ", nx,ny,nz);

   /* wave modeling space */
   int nxyz;
   nxyz=nx*ny*nz;

   c11.resize(nxyz);
   c22.resize(nxyz);
   c33.resize(nxyz);
   c12.resize(nxyz);
   c13.resize(nxyz);
   c23.resize(nxyz);
   c44.resize(nxyz);
   c55.resize(nxyz);
   c66.resize(nxyz);
   phx.resize(nxyz);
   phz.resize(nxyz);
   phy.resize(nxyz);
 
   c11r>>c11;
   c22r>>c22;
   c33r>>c33;
   c12r>>c12;
   c13r>>c13;
   c23r>>c23;
   c44r>>c44;
   c55r>>c55;
   c66r>>c66;
   phix>>phx;
   phiy>>phy;
   phiz>>phz;

  
   for(int i=0;i<nxyz;i++){
      phx[i] *= SF_PI/180.0;
      phy[i] *= SF_PI/180.0;
      phz[i] *= SF_PI/180.0;
   } 

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nx=%d ny=%d nz=%d nxyz=%d nk=%d",nx,ny,nz,nxyz,nk);

   float dkz,dkx,dky,kz0,kx0,ky0;

   dkx=2*SF_PI/dx/nx;
   dky=2*SF_PI/dy/ny;
   dkz=2*SF_PI/dz/nz;

   kx0=-SF_PI/dx;
   ky0=-SF_PI/dy;
   kz0=-SF_PI/dz;

   rkx.resize(nk);
   rky.resize(nk);
   rkz.resize(nk);
   rkk.resize(nk);

   double kx, ky, kz, k2, rk;
   int    i=0, j=0, k=0, ix, iy, iz;
   
   for(iy=0; iy < nky; iy++)
   {
     ky = ky0+iy*dky;
     if(ky==0.0) ky=0.0000000001*dky;

     for(ix=0; ix < nkx; ix++)
     {
       kx = kx0+ix*dkx;
       if(kx==0.0) kx=0.0000000001*dkx;

         for (iz=0; iz < nkz; iz++)
         {
            kz = kz0+iz*dkz;
            if(kz==0.0) kz=0.0000000001*dkz;

            k2 = ky*ky+kx*kx+kz*kz;
            rk = sqrt(k2);

            rky[i] = ky/rk;
            rkx[i] = kx/rk;
            rkz[i] = kz/rk;
            rkk[i] = rk;
            i++;
         }
      }
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for prereparing for low-rank decomp: %f(second)",timespent);

   /*****************************************************************************
   *  Calculating polarization deviation operator for wave-mode separation
   * ***************************************************************************/
   vector<int> md(nxyz), nd(nk);
   for (k=0; k < nxyz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /*****************************************************************************
                    low rank approximation for vector decomposition
   * ***************************************************************************/
   int   m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, m2xy, n2xy, m2xz, n2xz, m2yz, n2yz;

   /********* Ax2-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpxx,eps,npk,lid,rid,mid) );
   m2xx=mid.m();
   n2xx=mid.n();
   sf_warning("lowrank: m2xx=%d n2xx=%d",m2xx, n2xx);

   float *ldataxx, *fmidxx, *rdataxx;

   fmidxx  = sf_floatalloc(m2xx*n2xx);
   ldataxx = sf_floatalloc(nxyz*m2xx);
   rdataxx = sf_floatalloc(n2xx*nk);

   map2d1d(fmidxx, mid, m2xx, n2xx);

   iC ( samplebpxx(md,lid,mat) );
   map2d1d(ldataxx, mat, nxyz, m2xx);

   iC ( samplebpxx(rid,nd,mat) );
   map2d1d(rdataxx, mat, n2xx, nk);

   /********* Ay2-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpyy,eps,npk,lid,rid,mid) );
   m2yy=mid.m();
   n2yy=mid.n();
   sf_warning("lowrank: m2yy=%d n2yy=%d",m2yy, n2yy);

   float *ldatayy, *fmidyy, *rdatayy;

   fmidyy  = sf_floatalloc(m2yy*n2yy);
   ldatayy = sf_floatalloc(nxyz*m2yy);
   rdatayy = sf_floatalloc(n2yy*nk);

   map2d1d(fmidyy, mid, m2yy, n2yy);

   iC ( samplebpyy(md,lid,mat) );
   map2d1d(ldatayy, mat, nxyz, m2yy);

   iC ( samplebpyy(rid,nd,mat) );
   map2d1d(rdatayy, mat, n2yy, nk);

   /********* Az2-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpzz,eps,npk,lid,rid,mid) );
   m2zz=mid.m();
   n2zz=mid.n();
   sf_warning("lowrank: m2zz=%d n2zz=%d",m2zz, n2zz);

   float *ldatazz, *fmidzz, *rdatazz;

   fmidzz  = sf_floatalloc(m2zz*n2zz);
   ldatazz = sf_floatalloc(nxyz*m2zz);
   rdatazz = sf_floatalloc(n2zz*nk);

   map2d1d(fmidzz, mid, m2zz, n2zz);

   iC ( samplebpzz(md,lid,mat) );
   map2d1d(ldatazz, mat, nxyz, m2zz);

   iC ( samplebpzz(rid,nd,mat) );
   map2d1d(rdatazz, mat, n2zz, nk);

   /********* AxAy-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpxy,eps,npk,lid,rid,mid) );
   m2xy=mid.m();
   n2xy=mid.n();
   sf_warning("lowrank: m2xy=%d n2xy=%d",m2xy, n2xy);

   float *ldataxy, *fmidxy, *rdataxy;

   fmidxy  = sf_floatalloc(m2xy*n2xy);
   ldataxy = sf_floatalloc(nxyz*m2xy);
   rdataxy = sf_floatalloc(n2xy*nk);

   map2d1d(fmidxy, mid, m2xy, n2xy);

   iC ( samplebpxy(md,lid,mat) );
   map2d1d(ldataxy, mat, nxyz, m2xy);

   iC ( samplebpxy(rid,nd,mat) );
   map2d1d(rdataxy, mat, n2xy, nk);

   /********* AxAz-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpxz,eps,npk,lid,rid,mid) );
   m2xz=mid.m();
   n2xz=mid.n();
   sf_warning("lowrank: m2xz=%d n2xz=%d",m2xz, n2xz);

   float *ldataxz, *fmidxz, *rdataxz;

   fmidxz  = sf_floatalloc(m2xz*n2xz);
   ldataxz = sf_floatalloc(nxyz*m2xz);
   rdataxz = sf_floatalloc(n2xz*nk);

   map2d1d(fmidxz, mid, m2xz, n2xz);

   iC ( samplebpxz(md,lid,mat) );
   map2d1d(ldataxz, mat, nxyz, m2xz);

   iC ( samplebpxz(rid,nd,mat) );
   map2d1d(rdataxz, mat, n2xz, nk);

   /********* AyAz-operator **********/
   iC( ddlowrank(nxyz,nk,samplebpyz,eps,npk,lid,rid,mid) );
   m2yz=mid.m();
   n2yz=mid.n();
   sf_warning("lowrank: m2yz=%d n2yz=%d",m2yz, n2yz);

   float *ldatayz, *fmidyz, *rdatayz;

   fmidyz  = sf_floatalloc(m2yz*n2yz);
   ldatayz = sf_floatalloc(nxyz*m2yz);
   rdatayz = sf_floatalloc(n2yz*nk);

   map2d1d(fmidyz, mid, m2yz, n2yz);

   iC ( samplebpyz(md,lid,mat) );
   map2d1d(ldatayz, mat, nxyz, m2yz);

   iC ( samplebpyz(rid,nd,mat) );
   map2d1d(rdatayz, mat, n2yz, nk);


   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   //  wavelet parameter for source definition */

   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

    float *px, *py, *pz, *x;
    px=sf_floatalloc(nxyz);
    py=sf_floatalloc(nxyz);
    pz=sf_floatalloc(nxyz);
    x=sf_floatalloc(nxyz);

    int iflag;

    sf_floatread(px, nxyz, Fx);
    sf_floatread(py, nxyz, Fy);
    sf_floatread(pz, nxyz, Fz);

    sf_file Fpx, Fpy, Fpz;
    Fpx= sf_output("out");
    Fpy= sf_output("ElasticPy");
    Fpz= sf_output("ElasticPz");

    puthead3x(Fpx, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
    puthead3x(Fpy, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
    puthead3x(Fpz, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);

    // separate qP wave 
    sf_warning("vector decomposition of P-wave based on lowrank decomp.");
    decomplowrank3dp(ldataxx, rdataxx, fmidxx,
                    ldatayy, rdatayy, fmidyy,
                    ldatazz, rdatazz, fmidzz,
                    ldataxy, rdataxy, fmidxy,
                    ldataxz, rdataxz, fmidxz,
                    ldatayz, rdatayz, fmidyz,
                    px, py, pz, ijkx, ijky, ijkz,
                    nx, ny, nz, nxyz, nk,
                    m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, 
                    m2xz, n2xz, m2xy, n2xy, m2yz, n2yz);

	sf_floatwrite(px, nxyz, Fpx);
	sf_floatwrite(py, nxyz, Fpy);
	sf_floatwrite(pz, nxyz, Fpz);

    free(px);
    free(py);
    free(pz);

    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes decomposition.: %f(second)",timespent);

    free(ldataxx);
    free(ldatayy);
    free(ldatazz);
    free(ldataxy);
    free(ldataxz);
    free(ldatayz);

    free(rdataxx);
    free(rdatayy);
    free(rdatazz);
    free(rdataxy);
    free(rdataxz);
    free(rdatayz);

    free(fmidxx);
    free(fmidyy);
    free(fmidzz);
    free(fmidxy);
    free(fmidxz);
    free(fmidyz);

    sf_warning("-------sucessful ending --------");
    exit(0);
}

/* Apy2  */
static int samplebpyy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upy=-Chr[7];
            }

            res(a,b) = upy*upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apx2  */
static int samplebpxx(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
            }

            res(a,b) = upx*upx;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apz2  */
static int samplebpzz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upz=-Chr[8];
            }

            res(a,b) = upz*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApyApx  */
static int samplebpxy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
                upy=-Chr[7];
            }

            res(a,b) = upx*upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApxApz  */
static int samplebpxz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
                upz=-Chr[8];
            }

            res(a,b) = upx*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApyApz  */
static int samplebpyz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double b11, b22, b33, b12, b13, b23, b44, b55, b66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
	b11 = c11[i];
	b22 = c22[i];
	b33 = c33[i];
	b12 = c12[i];
	b13 = c13[i];
	b23 = c23[i];
	b44 = c44[i];
	b55 = c55[i];
	b66 = c66[i];

        double cosphiz=cos(phz[i]);
        double sinphiz=sin(phz[i]);
        double cosphiy=cos(phy[i]);
        double sinphiy=sin(phy[i]);
        double cosphix=cos(phx[i]);
        double sinphix=sin(phx[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cosphiy*cosphiz*kx0 + (cosphix*sinphiz + sinphix*sinphiy*cosphiz)*ky0 + (sinphix*sinphiz - cosphix*sinphiy*cosphiz)*kz0;
            ky = -cosphiy*sinphiz*kx0 + (cosphix*cosphiz - sinphix*sinphiy*sinphiz)*ky0 + (sinphix*cosphiz + cosphix*sinphiy*sinphiz)*kz0;
            kz = sinphiy*kx0 - sinphix*cosphiy*ky0 + cosphix*cosphiy*kz0;

            a11= b11*kx*kx + b66*ky*ky + b55*kz*kz;
            a22= b66*kx*kx + b22*ky*ky + b44*kz*kz;
            a33= b55*kx*kx + b44*ky*ky + b33*kz*kz;
            a12= (b12+b66)*kx*ky;
            a13= (b13+b55)*kx*kz;
            a23= (b23+b44)*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upy=-Chr[7];
                upz=-Chr[8];
            }

            res(a,b) = upy*upz;
              
         }// b loop
    }// a loop

    return 0;
}

static void map2d1d(float *d, DblNumMat mat, int m, int n)
{
   int i, j, k;
   k=0;
   for (i=0; i < m; i++)
   for (j=0; j < n; j++)
   {
        d[k] = (float)mat(i,j);
        k++;
   }

}

