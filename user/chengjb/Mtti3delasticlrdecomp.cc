/* 3-D three-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and vector decomposition using low-rank symbol approximation in 3D TTI media.

   Authors: Jiubing Cheng (Tongji University) and Sergey Fomel (The University of Texas at Austin)
     
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
#include <assert.h>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

/* prepared head files by myself */
#include "_cjb.h"
#include "_lapack.h"

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
#include "eigen3x3.h"
//#include "decomplowrank3d.h"
}

static std::valarray<float> vp, vs, ep, de, ga, th, ph;

static std::valarray<double> rkx, rky, rkz, rkk;

/* for low-rank decomp. */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int sampleyy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplexy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int sampleyz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplex2y2(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplex2z2(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int sampley2z2(vector<int>& rs, vector<int>& cs, DblNumMat& res);

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

   int   ns;
   float dt;

   par.get("ns",ns);
   par.get("dt",dt);

   int ireconstruct;   // flag for reconstruct the W matrix or not
   int iflagvti;       // 1: for VTI 
   par.get("iflagvti",iflagvti);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("iflagvti=%d ",iflagvti);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), gam("gam"), the("the"), phai("phai");

   /* Read/Write axes */
   int nxv, nyv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("n3",nyv);

   float az, ax, ay;
   vp0.get("o1",az);
   vp0.get("o2",ax);
   vp0.get("o3",ay);

   float fx, fy, fz;
   fx=ax*1000.0;
   fy=ay*1000.0;
   fz=az*1000.0;

   float dx, dy, dz;
   vp0.get("d1",az);
   vp0.get("d2",ax);
   vp0.get("d3",ay);
   dz = az*1000.0;
   dx = ax*1000.0;
   dy = ay*1000.0;

   /* wave modeling space */
   int nx, ny, nz, nxyz;
   nx=nxv;
   ny=nyv;
   nz=nzv;
   nxyz=nx*ny*nz;

   vp.resize(nxyz);
   vs.resize(nxyz);
   ep.resize(nxyz);
   de.resize(nxyz);
   ga.resize(nxyz);
   th.resize(nxyz);
   ph.resize(nxyz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;
   gam>>ga;
   the>>th;
   phai>>ph;
  
   for(int i=0;i<nxyz;i++){
      th[i] *= PI/180.0;
      ph[i] *= PI/180.0;
   } 

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nx=%d ny=%d nz=%d nxyz=%d nk=%d",nx,ny,nz,nxyz,nk);

   float dkz,dkx,dky,kz0,kx0,ky0;

   dkx=2*PI/dx/nx;
   dky=2*PI/dy/ny;
   dkz=2*PI/dz/nz;

   kx0=-PI/dx;
   ky0=-PI/dy;
   kz0=-PI/dz;

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
   int   m2x2y2, n2x2y2, m2x2z2, n2x2z2, m2y2z2, n2y2z2;

   /********* Ax2-operator **********/
   iC( ddlowrank(nxyz,nk,samplexx,eps,npk,lid,rid,mid) );
   m2xx=mid.m();
   n2xx=mid.n();
   sf_warning("lowrank: m2xx=%d n2xx=%d",m2xx, n2xx);

   float *ldataxx, *fmidxx, *rdataxx;

   fmidxx  = sf_floatalloc(m2xx*n2xx);
   ldataxx = sf_floatalloc(nxyz*m2xx);
   rdataxx = sf_floatalloc(n2xx*nk);

   map2d1d(fmidxx, mid, m2xx, n2xx);

   iC ( samplexx(md,lid,mat) );
   map2d1d(ldataxx, mat, nxyz, m2xx);

   iC ( samplexx(rid,nd,mat) );
   map2d1d(rdataxx, mat, n2xx, nk);

   /********* Ay2-operator **********/
   iC( ddlowrank(nxyz,nk,sampleyy,eps,npk,lid,rid,mid) );
   m2yy=mid.m();
   n2yy=mid.n();
   sf_warning("lowrank: m2yy=%d n2yy=%d",m2yy, n2yy);

   float *ldatayy, *fmidyy, *rdatayy;

   fmidyy  = sf_floatalloc(m2yy*n2yy);
   ldatayy = sf_floatalloc(nxyz*m2yy);
   rdatayy = sf_floatalloc(n2yy*nk);

   map2d1d(fmidyy, mid, m2yy, n2yy);

   iC ( sampleyy(md,lid,mat) );
   map2d1d(ldatayy, mat, nxyz, m2yy);

   iC ( sampleyy(rid,nd,mat) );
   map2d1d(rdatayy, mat, n2yy, nk);

   /********* Az2-operator **********/
   iC( ddlowrank(nxyz,nk,samplezz,eps,npk,lid,rid,mid) );
   m2zz=mid.m();
   n2zz=mid.n();
   sf_warning("lowrank: m2zz=%d n2zz=%d",m2zz, n2zz);

   float *ldatazz, *fmidzz, *rdatazz;

   fmidzz  = sf_floatalloc(m2zz*n2zz);
   ldatazz = sf_floatalloc(nxyz*m2zz);
   rdatazz = sf_floatalloc(n2zz*nk);

   map2d1d(fmidzz, mid, m2zz, n2zz);

   iC ( samplezz(md,lid,mat) );
   map2d1d(ldatazz, mat, nxyz, m2zz);

   iC ( samplezz(rid,nd,mat) );
   map2d1d(rdatazz, mat, n2zz, nk);

   /********* AxAy-operator **********/
   iC( ddlowrank(nxyz,nk,samplexy,eps,npk,lid,rid,mid) );
   m2xy=mid.m();
   n2xy=mid.n();
   sf_warning("lowrank: m2xy=%d n2xy=%d",m2xy, n2xy);

   float *ldataxy, *fmidxy, *rdataxy;

   fmidxy  = sf_floatalloc(m2xy*n2xy);
   ldataxy = sf_floatalloc(nxyz*m2xy);
   rdataxy = sf_floatalloc(n2xy*nk);

   map2d1d(fmidxy, mid, m2xy, n2xy);

   iC ( samplexy(md,lid,mat) );
   map2d1d(ldataxy, mat, nxyz, m2xy);

   iC ( samplexy(rid,nd,mat) );
   map2d1d(rdataxy, mat, n2xy, nk);

   /********* AxAz-operator **********/
   iC( ddlowrank(nxyz,nk,samplexz,eps,npk,lid,rid,mid) );
   m2xz=mid.m();
   n2xz=mid.n();
   sf_warning("lowrank: m2xz=%d n2xz=%d",m2xz, n2xz);

   float *ldataxz, *fmidxz, *rdataxz;

   fmidxz  = sf_floatalloc(m2xz*n2xz);
   ldataxz = sf_floatalloc(nxyz*m2xz);
   rdataxz = sf_floatalloc(n2xz*nk);

   map2d1d(fmidxz, mid, m2xz, n2xz);

   iC ( samplexz(md,lid,mat) );
   map2d1d(ldataxz, mat, nxyz, m2xz);

   iC ( samplexz(rid,nd,mat) );
   map2d1d(rdataxz, mat, n2xz, nk);

   /********* AyAz-operator **********/
   iC( ddlowrank(nxyz,nk,sampleyz,eps,npk,lid,rid,mid) );
   m2yz=mid.m();
   n2yz=mid.n();
   sf_warning("lowrank: m2yz=%d n2yz=%d",m2yz, n2yz);

   float *ldatayz, *fmidyz, *rdatayz;

   fmidyz  = sf_floatalloc(m2yz*n2yz);
   ldatayz = sf_floatalloc(nxyz*m2yz);
   rdatayz = sf_floatalloc(n2yz*nk);

   map2d1d(fmidyz, mid, m2yz, n2yz);

   iC ( sampleyz(md,lid,mat) );
   map2d1d(ldatayz, mat, nxyz, m2yz);

   iC ( sampleyz(rid,nd,mat) );
   map2d1d(rdatayz, mat, n2yz, nk);


   /********* Ax2+Ay2-operator **********/
   iC( ddlowrank(nxyz,nk,samplex2y2,eps,npk,lid,rid,mid) );
   m2x2y2=mid.m();
   n2x2y2=mid.n();
   sf_warning("lowrank: m2x2y2=%d n2x2y2=%d",m2x2y2, n2x2y2);

   float *ldatax2y2, *fmidx2y2, *rdatax2y2;

   fmidx2y2  = sf_floatalloc(m2x2y2*n2x2y2);
   ldatax2y2 = sf_floatalloc(nxyz*m2x2y2);
   rdatax2y2 = sf_floatalloc(n2x2y2*nk);

   map2d1d(fmidx2y2, mid, m2x2y2, n2x2y2);

   iC ( samplex2y2(md,lid,mat) );
   map2d1d(ldatax2y2, mat, nxyz, m2x2y2);

   iC ( samplex2y2(rid,nd,mat) );
   map2d1d(rdatax2y2, mat, n2x2y2, nk);

   /********* Ax2+Az2-operator **********/
   iC( ddlowrank(nxyz,nk,samplex2z2,eps,npk,lid,rid,mid) );
   m2x2z2=mid.m();
   n2x2z2=mid.n();
   sf_warning("lowrank: m2x2z2=%d n2x2z2=%d",m2x2z2, n2x2z2);

   float *ldatax2z2, *fmidx2z2, *rdatax2z2;

   fmidx2z2  = sf_floatalloc(m2x2z2*n2x2z2);
   ldatax2z2 = sf_floatalloc(nxyz*m2x2z2);
   rdatax2z2 = sf_floatalloc(n2x2z2*nk);

   map2d1d(fmidx2z2, mid, m2x2z2, n2x2z2);

   iC ( samplex2z2(md,lid,mat) );
   map2d1d(ldatax2z2, mat, nxyz, m2x2z2);

   iC ( samplex2z2(rid,nd,mat) );
   map2d1d(rdatax2z2, mat, n2x2z2, nk);

   /********* Ay2+Az2-operator **********/
   iC( ddlowrank(nxyz,nk,sampley2z2,eps,npk,lid,rid,mid) );
   m2y2z2=mid.m();
   n2y2z2=mid.n();
   sf_warning("lowrank: m2y2z2=%d n2y2z2=%d",m2y2z2, n2y2z2);

   float *ldatay2z2, *fmidy2z2, *rdatay2z2;

   fmidy2z2  = sf_floatalloc(m2y2z2*n2y2z2);
   ldatay2z2 = sf_floatalloc(nxyz*m2y2z2);
   rdatay2z2 = sf_floatalloc(n2y2z2*nk);

   map2d1d(fmidy2z2, mid, m2y2z2, n2y2z2);

   iC ( sampley2z2(md,lid,mat) );
   map2d1d(ldatay2z2, mat, nxyz, m2y2z2);

   iC ( sampley2z2(rid,nd,mat) );
   map2d1d(rdatay2z2, mat, n2y2z2, nk);

   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   //  wavelet parameter for source definition */
   /*

   int nxpad, nypad, nzpad;
   nxpad=nx+2*M;
   nypad=ny+2*M;
   nzpad=nz+2*M;

   sf_warning("fx=%f fy=%f fz=%f dx=%f dy=%f dz=%f",fx,fy,fz,dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d nxpad=%d nypad=%d nzpad=%d", nx,ny,nz,nxpad,nypad,nzpad);

   int mm=2*M+1;

   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

   // source definition 
   int ixs, iys, izs, ixms, iyms, izms;

   // setup I/O files 
   iRSF Elasticx("Elasticx"),Elasticy("Elasticy"),Elasticz("Elasticz");

   Elasticx.get("n1",nz);
   Elasticx.get("n2",nx);
   Elasticx.get("n3",ny);
   Elasticx.get("d1",dz/1000);
   Elasticx.get("d2",dx/1000);
   Elasticx.get("d3",dy/1000);
   Elasticx.get("o1",fz/1000);
   Elasticx.get("o2",fx/1000);
   Elasticx.get("o3",fy/1000);

    float **p=sf_floatalloc3(nz, nx, ny);
    float **q=sf_floatalloc3(nz, nx, ny);
    float **r=sf_floatalloc3(nz, nx, ny);

    zero3float(p, nz, nx, ny);
    zero3float(q, nz, nx, ny);
    zero3float(r, nz, nx, ny);

    std::valarray<float> x(nz);

    float *pp, *qq, *rr;
    pp=sf_floatalloc(nxyz);
    qq=sf_floatalloc(nxyz);
    rr=sf_floatalloc(nxyz);

    int ii, jj, im, jm;

    int iflag=0;
    k=0;
    for(l=0;l<ny;l++)
    {
       lm=l+M;
       for(i=0;i<nx;i++)
       {
         im=i+M;
         for(j=0;j<nz;j++)
	 {
            jm=j+M;
            pp[k] = p[im][jm];
            qq[k] = q[im][jm];
            rr[k] = r[im][jm];

            k++;      
         }
      }
     }// i loop

    k=0;
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++){
       for(iz=0;iz<nz;iz++){
          x[iz] = pp[k] + qq[k] + rr[k];
          k++;
       }
       sf_floatwrite(x,nz,Fo1);
    }

    // separate qP wave 
    sf_warning("vector decomposition of P-wave based on lowrank decomp.");
    decomplowrank3dp(ldataxx, rdataxx, fmidxx,
                    ldatayy, rdatayy, fmidyy,
                    ldatazz, rdatazz, fmidzz,
                    ldataxy, rdataxy, fmidxy,
                    ldataxz, rdataxz, fmidxz,
                    ldatayz, rdatayz, fmidyz,
                    p3, q3, r3, px, py, pz, ijkx, ijky, ijkz,
                    nx, ny, nz, nxyz, nk, M,
                    m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, 
                    m2xz, n2xz, m2xy, n2xy, m2yz, n2yz);

    for(i=0;i<nxyz;i++)
         x[i] = px[i];
    ElasticPx<<x;

    for(i=0;i<nxyz;i++)
         x[i] = py[i];
    ElasticPy<<x;

    for(i=0;i<nxyz;i++)
         x[i] = pz[i];
    ElasticPz<<x;

    // separate qS wave (SV + SH) 
    sf_warning("vector decomposition of S-wave based on lowrank decomp.");
    decomplowrank3ds(ldatax2y2, rdatax2y2, fmidx2y2,
                    ldatax2z2, rdatax2z2, fmidx2z2,
                    ldatay2z2, rdatay2z2, fmidy2z2,
                    ldataxy, rdataxy, fmidxy,
                    ldataxz, rdataxz, fmidxz,
                    ldatayz, rdatayz, fmidyz,
                    p3, q3, r3, px, py, pz, ijkx, ijky, ijkz,
                    nx, ny, nz, nxyz, nk, M,
                    m2x2y2, n2x2y2, m2x2z2, n2x2z2, m2y2z2, n2y2z2,
                    m2xz, n2xz, m2xy, n2xy, m2yz, n2yz);

    for(i=0;i<nxyz;i++)
         x[i] = px[i];
    ElasticSx<<x;

    for(i=0;i<nxyz;i++)
         x[i] = py[i];
    ElasticSy<<x;

    for(i=0;i<nxyz;i++)
         x[i] = pz[i];
    ElasticSz<<x;
    }

    free(pp);
    free(qq);
    free(rr);

    free(**p);
    free(**q);
    free(**r);

    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);
    */

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

    free(fmidx2y2);
    free(fmidx2z2);
    free(fmidy2z2);

    sf_warning("-------sucessful ending --------");
    exit(0);
}

/* Apy2  */
static int sampleyy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upy*upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apx2  */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upx*upx;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apz2  */
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upz*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApyApx  */
static int samplexy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upx*upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApxApz  */
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upx*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApyApz  */
static int sampleyz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upy*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apy2+Apx2  */
static int samplex2y2(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upx*upx+ upy*upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apx2 + Apz2  */
static int samplex2z2(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upx*upx +upz*upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* Apy2 + Apz2 */
static int sampley2z2(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

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

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_ph*kx0 - sin_ph*ky0 + sin_th*cos_ph*kz0;
            ky = cos_th*sin_ph*kx0 + cos_ph*ky0 + sin_th*sin_ph*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

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
                upz=-Chr[8];
            }

            res(a,b) = upy*upy + upz*upz;
              
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

