/* 2-D two-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and vector decomposition based on lowrank approximation in TTI media.

   Authors: Jiubing Cheng
     
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

/*#include "_fd.h" 
 * when _fd.h has declared "m", we met compiling error at "m2=mid.m()";
 * "expected unqualified-id before numeric constant" 
*/
#ifndef M
#define M 5   /* 10th-order finite-difference: accuracy is 2*M */
#endif
#ifndef Mix
#define Mix 5   /* order of finite-difference for mix-derivative (2*mix) */
#endif

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
#include "fdcoef.h"
#include "eigen2x2.h"
#include "fwpttielastic.h"
#include "decomplowrank.h"
#include "seplowrank.h"  // for sub-program: reconstruct(...)
}

static std::valarray<float> vp, vs, ep, de, th;

static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation and wave vector decomposition operator based on low-rank decomp. */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static void polxp2dtti(float **ap, int nx, int nz, int im);
static void polzp2dtti(float **ap, int nx, int nz, int im);

static void map2d1d(float *d, DblNumMat mat, int m, int n);
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
   float dt=0.004, xrec=1000.0, zrec=1000.0;

   par.get("ns",ns);
   par.get("dt",dt);

  int ireconstruct;   // flag for reconstruct the W matrix or not
  par.get("ireconstruct",ireconstruct);
   par.get("xrec",xrec);
   par.get("zrec",zrec);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("ireconstruct=%d ",ireconstruct);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), the("the");

   /* Read/Write axes */
   int nxv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);

   float az, ax;
   vp0.get("o1",az);
   vp0.get("o2",ax);

   float fx, fz;
   fx=ax*1000.0;
   fz=az*1000.0;

   float dx, dz;
   vp0.get("d1",az);
   vp0.get("d2",ax);
   dz = az*1000.0;
   dx = ax*1000.0;

   /* location of operator check */
   int ixrec=(xrec-fx)/dx;
   int izrec=(zrec-fz)/dz;
   sf_warning("xrec=%f ",xrec);
   sf_warning("zrec=%f ",zrec);
   sf_warning("ixrec=%d ",ixrec);
   sf_warning("izrec=%d ",izrec);

   /* wave modeling space */
   int nx, nz, nxz;
   nx=nxv;
   nz=nzv;
   nxz=nx*nz;

   vp.resize(nxz);
   vs.resize(nxz);
   ep.resize(nxz);
   de.resize(nxz);
   th.resize(nxz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;
   the>>th;

   for(int i=0;i<nxz;i++)
      th[i] *= SF_PI/180.0;

   /* Fourier spectra demension */
   int nkz,nkx,nk;
   nkx=nx;
   nkz=nz;
   nk = nkx*nkz;

   float dkz,dkx,kz0,kx0;

   dkx=2*SF_PI/dx/nx;
   dkz=2*SF_PI/dz/nz;

   kx0=-SF_PI/dx;
   kz0=-SF_PI/dz;

   sinx.resize(nk);
   cosx.resize(nk);
   rkk.resize(nk);

   double kx, kz, k2, rk;
   int    i=0, j=0, k=0, ix, iz;
   
   for(ix=0; ix < nkx; ix++)
   {
       kx = kx0+ix*dkx;
       if(kx==0.0) kx=0.0000000001*dkx;

       for (iz=0; iz < nkz; iz++)
       {
            kz = kz0+iz*dkz;
            if(kz==0.0) kz=0.0000000001*dkz;

            k2 = kx*kx+kz*kz;
            rk = sqrt(k2);

            sinx[i] = kx/rk;
            cosx[i] = kz/rk;
            rkk[i] = rk;

            i++;
       }
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for prereparing for low-rank decomp: %f(second)",timespent);

   /*****************************************************************************
   *  Calculating vector decomposition operator for P-wave
   * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   int   m2xx, n2xx, m2xz, n2xz, m2zz, n2zz;
   float *ldataxx, *fmidxx, *rdataxx;
   float *ldataxz, *fmidxz, *rdataxz;
   float *ldatazz, *fmidzz, *rdatazz;

   /** Apx2 ****/
   iC( ddlowrank(nxz,nk,samplexx,eps,npk,lid,rid,mid) );
   m2xx=mid.m();
   n2xx=mid.n();
   sf_warning("m2xx=%d n2xx=%d",m2xx, n2xx);

   fmidxx  = sf_floatalloc(m2xx*n2xx);
   ldataxx = sf_floatalloc(nxz*m2xx);
   rdataxx = sf_floatalloc(n2xx*nk);

   map2d1d(fmidxx, mid, m2xx, n2xx);

   iC ( samplexx(md,lid,mat) );
   map2d1d(ldataxx, mat, nxz, m2xx);

   iC ( samplexx(rid,nd,mat) );
   map2d1d(rdataxx, mat, n2xx, nk);

   /********* reconsturct W-matrix for checking **********/
   int im=ixrec*nz+izrec;
   float *w;
   float **apx, **apz;
   std::valarray<float> w1(nk);
   if(ireconstruct==1)
   {
      w = sf_floatalloc(nk);

      sf_warning("reconstruct x-component operators based on low-rank approx.");
      reconstruct1(w, ldataxx, fmidxx, rdataxx, nxz, nk, m2xx, n2xx, im);

      /* error check between accurate and low-rank approx. */
      apx=sf_floatalloc2(nz, nx);

      oRSF Errdecxp("Errdecxp");

      Errdecxp.put("n1",nkz);
      Errdecxp.put("n2",nkx);
      Errdecxp.put("d1",dkz);
      Errdecxp.put("d2",dkx);
      Errdecxp.put("o1",kz0);
      Errdecxp.put("o2",kx0);

      Errdecxp.put("label1","kz");
      Errdecxp.put("label2","kx");
      Errdecxp.put("unit1","2*pi/m");
      Errdecxp.put("unit2","2*pi/m");

      //calculate P-wave's polarization operators directly for X-component
      polxp2dtti(apx, nx, nz, im);

      float sum=0.0;
      int k=0;
      for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
      {
           float err=w[k]-apx[i][j]*apx[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
       }
       sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);

       Errdecxp << w1;

      oRSF Decompxp("Decompxp");

      Decompxp.put("n1",nkz);
      Decompxp.put("n2",nkx);
      Decompxp.put("d1",dkz);
      Decompxp.put("d2",dkx);
      Decompxp.put("o1",kz0);
      Decompxp.put("o2",kx0);

      Decompxp.put("label1","kz");
      Decompxp.put("label2","kx");
      Decompxp.put("unit1","2*pi/m");
      Decompxp.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++)
         w1[i]=w[i];

      Decompxp << w1;

   }//ireconstruct==1

   /** ApxApz ****/
   iC( ddlowrank(nxz,nk,samplexz,eps,npk,lid,rid,mid) );
   m2xz=mid.m();
   n2xz=mid.n();
   sf_warning("m2xz=%d n2xz=%d",m2xz, n2xz);

   fmidxz  = sf_floatalloc(m2xz*n2xz);
   ldataxz = sf_floatalloc(nxz*m2xz);
   rdataxz = sf_floatalloc(n2xz*nk);

   map2d1d(fmidxz, mid, m2xz, n2xz);

   iC ( samplexz(md,lid,mat) );
   map2d1d(ldataxz, mat, nxz, m2xz);

   iC ( samplexz(rid,nd,mat) );
   map2d1d(rdataxz, mat, n2xz, nk);

   /********* reconsturct W-matrix for checking **********/
   if(ireconstruct==1)
   {
      sf_warning("reconstruct xz-component operators based on low-rank approx.");
      reconstruct1(w, ldataxz, fmidxz, rdataxz, nxz, nk, m2xz, n2xz, im);

      /* error check between accurate and low-rank approx. */
      apz=sf_floatalloc2(nz, nx);

      oRSF Errdecxzp("Errdecxzp");

      Errdecxzp.put("n1",nkz);
      Errdecxzp.put("n2",nkx);
      Errdecxzp.put("d1",dkz);
      Errdecxzp.put("d2",dkx);
      Errdecxzp.put("o1",kz0);
      Errdecxzp.put("o2",kx0);

      Errdecxzp.put("label1","kz");
      Errdecxzp.put("label2","kx");
      Errdecxzp.put("unit1","2*pi/m");
      Errdecxzp.put("unit2","2*pi/m");

      //calculate P-wave's polarization operators directly for X- and Z-component
      polxp2dtti(apx, nx, nz, im);
      polzp2dtti(apz, nx, nz, im);

      float sum=0.0;
      int k=0;
      for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
      {
           float err=w[k]-apx[i][j]*apz[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
      }
      sf_warning("XZ-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);

      Errdecxzp << w1;

      oRSF Decompxzp("Decompxzp");

      Decompxzp.put("n1",nkz);
      Decompxzp.put("n2",nkx);
      Decompxzp.put("d1",dkz);
      Decompxzp.put("d2",dkx);
      Decompxzp.put("o1",kz0);
      Decompxzp.put("o2",kx0);

      Decompxzp.put("label1","kz");
      Decompxzp.put("label2","kx");
      Decompxzp.put("unit1","2*pi/m");
      Decompxzp.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++)
         w1[i]=w[i];

      Decompxzp << w1;

      free(*apx);
   }//ireconstruct==1

   /** Apz2 ****/
   iC( ddlowrank(nxz,nk,samplezz,eps,npk,lid,rid,mid) );
   m2zz=mid.m();
   n2zz=mid.n();
   sf_warning("m2zz=%d n2zz=%d",m2zz, n2zz);

   fmidzz  = sf_floatalloc(m2zz*n2zz);
   ldatazz = sf_floatalloc(nxz*m2zz);
   rdatazz = sf_floatalloc(n2zz*nk);

   map2d1d(fmidzz, mid, m2zz, n2zz);

   iC ( samplezz(md,lid,mat) );
   map2d1d(ldatazz, mat, nxz, m2zz);

   iC ( samplezz(rid,nd,mat) );
   map2d1d(rdatazz, mat, n2zz, nk);

   if(ireconstruct==1)
   {
      sf_warning("reconstruct z-component operators based on low-rank approx.");
      reconstruct1(w, ldatazz, fmidzz, rdatazz, nxz, nk, m2zz, n2zz, im);

      /* error check between accurate and low-rank approx. */
      oRSF Errdeczp("Errdeczp");

      Errdeczp.put("n1",nkz);
      Errdeczp.put("n2",nkx);
      Errdeczp.put("d1",dkz);
      Errdeczp.put("d2",dkx);
      Errdeczp.put("o1",kz0);
      Errdeczp.put("o2",kx0);

      Errdeczp.put("label1","kz");
      Errdeczp.put("label2","kx");
      Errdeczp.put("unit1","2*pi/m");
      Errdeczp.put("unit2","2*pi/m");

      //calculate P-wave's polarization operators directly for Z-component
      polzp2dtti(apz, nx, nz, im);

      float sum=0.0;
      int k=0;
      for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
      {
           float err=w[k]-apz[i][j]*apz[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
       }
       sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);

       Errdeczp << w1;

      oRSF Decompzp("Decompzp");

      Decompzp.put("n1",nkz);
      Decompzp.put("n2",nkx);
      Decompzp.put("d1",dkz);
      Decompzp.put("d2",dkx);
      Decompzp.put("o1",kz0);
      Decompzp.put("o2",kx0);

      Decompzp.put("label1","kz");
      Decompzp.put("label2","kx");
      Decompzp.put("unit1","2*pi/m");
      Decompzp.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++)
         w1[i]=w[i];

      Decompzp << w1;

      free(w);
      free(*apz);
   }//ireconstruct==1

   /****************End of Calculating vector Decomposition Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   /*  wavelet parameter for source definition */
   float A, f0, t0;
   f0=30.0;                  
   t0=0.04;                  
   A=1.0;

   int nxpad, nzpad;
   nxpad=nx+2*M;
   nzpad=nz+2*M;

   sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);
   sf_warning("nx=%d nz=%d nxpad=%d nzpad=%d", nx,nz,nxpad,nzpad);

   int mm=2*M+1;

   float dt2=dt*dt;

   /* source definition */
   int ixs, izs, ixms, izms;
   ixs=nxv/2;
   izs=nzv/2;
   ixms=ixs+M;  /* source's x location */
   izms=izs+M;  /* source's z-location */

   float xs, zs;
   xs=fx+ixs*dx;
   zs=fz+izs*dz;
   sf_warning("source location: (%f, %f)",xs, zs);

   /* setup I/O files */
   oRSF ElasticX("out"),ElasticZ("ElasticZ");
   oRSF ElasticPx("ElasticPx");
   oRSF ElasticPz("ElasticPz");
   oRSF ElasticSVx("ElasticSVx");
   oRSF ElasticSVz("ElasticSVz");

   ElasticX.put("n1",nkz);
   ElasticX.put("n2",nkx);
   ElasticX.put("d1",dz/1000);
   ElasticX.put("d2",dx/1000);
   ElasticX.put("o1",fz/1000);
   ElasticX.put("o2",fx/1000);

   ElasticZ.put("n1",nkz);
   ElasticZ.put("n2",nkx);
   ElasticZ.put("d1",dz/1000);
   ElasticZ.put("d2",dx/1000);
   ElasticZ.put("o1",fz/1000);
   ElasticZ.put("o2",fx/1000);

   ElasticPx.put("n1",nkz);
   ElasticPx.put("n2",nkx);
   ElasticPx.put("d1",dz/1000);
   ElasticPx.put("d2",dx/1000);
   ElasticPx.put("o1",fz/1000);
   ElasticPx.put("o2",fx/1000);

   ElasticPz.put("n1",nkz);
   ElasticPz.put("n2",nkx);
   ElasticPz.put("d1",dz/1000);
   ElasticPz.put("d2",dx/1000);
   ElasticPz.put("o1",fz/1000);
   ElasticPz.put("o2",fx/1000);

   ElasticSVx.put("n1",nkz);
   ElasticSVx.put("n2",nkx);
   ElasticSVx.put("d1",dz/1000);
   ElasticSVx.put("d2",dx/1000);
   ElasticSVx.put("o1",fz/1000);
   ElasticSVx.put("o2",fx/1000);

   ElasticSVz.put("n1",nkz);
   ElasticSVz.put("n2",nkx);
   ElasticSVz.put("d1",dz/1000);
   ElasticSVz.put("d2",dx/1000);
   ElasticSVz.put("o1",fz/1000);
   ElasticSVz.put("o2",fx/1000);

    float *coeff_2dx=sf_floatalloc(mm);
    float *coeff_2dz=sf_floatalloc(mm);
    float *coeff_1dx=sf_floatalloc(mm);
    float *coeff_1dz=sf_floatalloc(mm);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dz,dz);
    coeff1dmix(coeff_1dx,dx);
    coeff1dmix(coeff_1dz,dz);

    float **p1=sf_floatalloc2(nzpad, nxpad);
    float **p2=sf_floatalloc2(nzpad, nxpad);
    float **p3=sf_floatalloc2(nzpad, nxpad);

    float **q1=sf_floatalloc2(nzpad, nxpad);
    float **q2=sf_floatalloc2(nzpad, nxpad);
    float **q3=sf_floatalloc2(nzpad, nxpad);

    zero2float(p1, nzpad, nxpad);
    zero2float(p2, nzpad, nxpad);
    zero2float(p3, nzpad, nxpad);

    zero2float(q1, nzpad, nxpad);
    zero2float(q2, nzpad, nxpad);
    zero2float(q3, nzpad, nxpad);

    sf_warning("==================================================");
    sf_warning("==  Porpagation Using Elastic anisotropic Eq.   ==");
    sf_warning("==================================================");

    float **vpp, **vss, **epp, **dee, **thee;

    vpp=sf_floatalloc2(nz,nx);
    vss=sf_floatalloc2(nz,nx);
    epp=sf_floatalloc2(nz,nx);
    dee=sf_floatalloc2(nz,nx);
    thee=sf_floatalloc2(nz,nx);

    k=0;
    for(i=0;i<nx;i++)
    for(j=0;j<nz;j++)
    {
       vpp[i][j]=vp[k];
       vss[i][j]=vs[k];
       epp[i][j]=ep[k];
       dee[i][j]=de[k];
       thee[i][j]=th[k];
       k++;
    }

    t4=clock();
    timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
    sf_warning("CPU time for preparing for wavefield modeling: %f(second)",timespent);

    std::valarray<float> x(nxz);
    float *px, *pz;
    px=sf_floatalloc(nxz);
    pz=sf_floatalloc(nxz);

    int ii, jj;

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

    for(int it=0;it<ns;it++)
    {
	float t=it*dt;

	if(it%50==0)
		sf_warning("Elastic: it= %d",it);

        // 2D exploding force source (e.g., Wu's PhD
        for(i=-1;i<=1;i++)
        for(j=-1;j<=1;j++)
        {
             if(fabs(i)+fabs(j)==2)
             {
                  p2[ixms+i][izms+j]+=i*Ricker(t, f0, t0, A);
                  q2[ixms+i][izms+j]+=j*Ricker(t, f0, t0, A);
             }
        }
        // 2D equil-energy force source (e.g., Wu's PhD)
        /*
        for(i=-1;i<=1;i++)
        for(j=-1;j<=1;j++)
        {
             if(fabs(i)+fabs(j)==2)
             {
                  if(i==-1&&j==1)  
                    q2[ixms+i][izms+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==-1&&j==-1) 
                   p2[ixms+i][izms+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==1&&j==1)  
                   p2[ixms+i][izms+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==1&&j==-1) 
                    q2[ixms+i][izms+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
             }
        }
        */

        /* fwpttielastic: forward-propagating using original elastic equation of displacement in TTI media*/
        fwpttielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
                      dx, dz, nx, nz, nxpad, nzpad, vpp, vss, epp, dee, thee);

        /******* output wavefields: component and divergence *******/
        if(it==ns-1)
	{
              k=0;
              for(i=0;i<nx;i++){
                int im=i+M; 
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   x[k] = px[k] = p3[im][jm];
                   k++;
                }
              }
              ElasticX<<x;

              k=0;
              for(i=0;i<nx;i++){
                int im=i+M; 
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   x[k] = pz[k] = q3[im][jm];
                   k++;
                }
              }
              ElasticZ<<x;

              t44=clock();
              timespent=(float)(t44-t4)/CLOCKS_PER_SEC;
              sf_warning("CPU time for wavefield modeling: %f(second)",timespent);

              /* separate qP wave  */
              sf_warning("vector decomposition of P-wave based on lowrank decomp."); 
              decomplowrank2d(ldataxx, rdataxx, fmidxx, 
                              ldataxz, rdataxz, fmidxz,                              
                              ldatazz, rdatazz, fmidzz,
                              px, pz, ijkx, ijkz,
                              nx, nz, nxz, nk, M, m2xx, n2xx, m2xz, n2xz, m2zz, n2zz);

              for(i=0;i<nxz;i++)
                 x[i] = px[i];
              ElasticPx<<x;

              for(i=0;i<nxz;i++)
                 x[i] = pz[i];
              ElasticPz<<x;
          
              k=0;
              for(i=0;i<nx;i++){
                int im=i+M;
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   px[k] = p3[im][jm];
                   pz[k] = q3[im][jm];
                   k++;
                }
              }
              /* separate qSV wave  */
              sf_warning("vector decomposition of SV-wave based on lowrank decomp.");
              decomplowrank2ds(ldataxx, rdataxx, fmidxx,
                               ldataxz, rdataxz, fmidxz,
                               ldatazz, rdatazz, fmidzz,
                               px, pz, ijkx, ijkz,
                               nx, nz, nxz, nk, M, m2xx, n2xx, m2xz, n2xz, m2zz, n2zz);

              for(i=0;i<nxz;i++)
                 x[i] = px[i];
              ElasticSVx<<x;

              for(i=0;i<nxz;i++)
                 x[i] = pz[i];
              ElasticSVz<<x;

              t5=clock();
              timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
              sf_warning("CPU time for vector decomposition.: %f(second)",timespent);

         }/* (it+1)%ntstep==0 */

         /**************************************/
 	 for(i=0,ii=M;i<nx;i++,ii++)
	    for(j=0,jj=M;j<nz;j++,jj++)
	    {
		p1[ii][jj]=p2[ii][jj];	
		p2[ii][jj]=p3[ii][jj];	

		q1[ii][jj]=q2[ii][jj];	
		q2[ii][jj]=q3[ii][jj];	
	    }

    }/* it loop */

    free(ldataxx);
    free(ldatazz);
    free(ldataxz);

    free(rdataxx);
    free(rdataxz);
    free(rdatazz);

    free(fmidxx);
    free(fmidxz);
    free(fmidzz);

    free(px);
    free(pz);

    free(*p1);
    free(*p2);
    free(*p3);
    free(*q1);
    free(*q2);
    free(*q3);

    free(*vpp);
    free(*vss);
    free(*epp);
    free(*dee);
    free(*thee);

    free(coeff_2dx);
    free(coeff_2dz);
    free(coeff_1dx);
    free(coeff_1dz);

    free(ijkx);
    free(ijkz);
exit(0);
}
/* Apx2 */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    double sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            if(s==0&&c==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            // rotatiing according to tilted symmetry axis
            sx=s*coss+c*sins;
            cx=c*coss-s*sins;
              
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u1*u1;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApxApz*/
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    double sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];

            if(s==0&&c==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            // rotatiing according to tilted symmetry axis
            sx=s*coss+c*sins;
            cx=c*coss-s*sins;
              
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u1*u2;
              
         }// b loop
    }// a loop

    return 0;
}
/* Apz2 */ 
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    double sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            if(s==0&&c==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            // rotatiing according to tilted symmetry axis
            sx=s*coss+c*sins;
            cx=c*coss-s*sins;
              
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u2*u2;
              
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

static void polxp2dtti(float **ap, int nx, int nz, int im)
{
        int    i, j, k;
        double sx, cx;

        double ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/
        double a[2][2];

        double vp2 = vp[im]*vp[im];
        double vs2 = vs[im]*vs[im];
        double ep2 = 1.0+2*ep[im];
        double de2 = 1.0+2*de[im];

        double u1, u2, c11, c33, c44, c13c44, a11, a12, a22;

        double coss=cos(th[im]);
        double sins=sin(th[im]);

        c33=vp2;
        c44=vs2;
        c11=ep2*c33;
        c13c44=sqrt((de2*c33-c44)*(c33-c44));

        for( i=0; i<nx ; i++ )
           for( j=0; j<nz ; j++)
                  ap[i][j]=0.0;

        k=0;
        for( i=0; i<nx ; i++ )
           for( j=0; j<nz ; j++)
           {
                sx=sinx[k];
                cx=cosx[k];
                k++;

                // rotatiing according to tilted symmetry axis
                double s=sx*coss+cx*sins;
                double c=cx*coss-sx*sins;

                a11= c11*s*s + c44*c*c;
                a12= c13c44*s*c;
                a22= c44*s*s + c33*c*c;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[0][0];
               u2=ve[0][1];
               /* get the closest direction to k */
               if(u1*s + u2*c <0) {
                  u1 = -u1;
               }
	       ap[i][j]=(float)(u1);
          } /* j loop */
}

static void polzp2dtti(float **ap, int nx, int nz, int im)
{
        int    i, j, k;
        double sx, cx;

        double ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/
        double a[2][2];

        //double f=1.0-vs2/vp2;

        double vp2 = vp[im]*vp[im];
        double vs2 = vs[im]*vs[im];
        double ep2 = 1.0+2*ep[im];
        double de2 = 1.0+2*de[im];

        double coss=cos(th[im]);
        double sins=sin(th[im]);

        double u1, u2, c11, c33, c44, c13c44, a11, a12, a22;

        c33=vp2;
        c44=vs2;
        c11=ep2*c33;
        c13c44=sqrt((de2*c33-c44)*(c33-c44));

        for( i=0; i<nx ; i++ )
           for( j=0; j<nz ; j++)
                  ap[i][j]=0.0;

        k=0;
        for( i=0; i<nx ; i++ )
           for( j=0; j<nz ; j++)
           {
                sx=sinx[k];
                cx=cosx[k];
                k++;

                // rotatiing according to tilted symmetry axis
                double s=sx*coss+cx*sins;
                double c=cx*coss-sx*sins;

                a11= c11*s*s + c44*c*c;
                a12= c13c44*s*c;
                a22= c44*s*s + c33*c*c;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[0][0];
               u2=ve[0][1];
               /* get the closest direction to k */
               if(u1*s + u2*c <0) {
                  u2 = -u2;
               }
	       ap[i][j]=(float)(u2);
          } /* j loop */
}

