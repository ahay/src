/* 2-D two-components wavefield modeling using original elastic anisotropic displacement 
  wave equation in VTI media.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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
#include "fwpvtielastic.h"
#include "seplowrank.h"
}

static std::valarray<float> vp, vs, ep, de;

static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplexs(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezs(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

static void polxp2dvti(float **ap, int nx, int nz, int im);
static void polzp2dvti(float **ap, int nx, int nz, int im);
static void polxs2dvti(float **ap, int nx, int nz, int im);
static void polzs2dvti(float **ap, int nx, int nz, int im);

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
   par.get("ireconstruct",ireconstruct);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("ireconstruct=%d ",ireconstruct);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del");

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

   /* wave modeling space */
   int nx, nz, nxz;
   nx=nxv;
   nz=nzv;
   nxz=nx*nz;

   vp.resize(nxz);
   vs.resize(nxz);
   ep.resize(nxz);
   de.resize(nxz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;

   /* Fourier spectra demension */
   int nkz,nkx,nk;
   nkx=nx;
   nkz=nz;
   nk = nkx*nkz;

   float dkz,dkx,kz0,kx0;

   dkx=2*PI/dx/nx;
   dkz=2*PI/dz/nz;

   kx0=-PI/dx;
   kz0=-PI/dz;

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
    *  Calculating polarization deviation operator for wave-mode separation
    * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /********* low rank decomposition p-wave, x-component **********/
   int   m2xp, n2xp, m2zp, n2zp;
   float *ldataxp, *fmidxp, *rdataxp;
   float *ldatazp, *fmidzp, *rdatazp;

   iC( ddlowrank(nxz,nk,samplexp,eps,npk,lid,rid,mid) );
   m2xp=mid.m();
   n2xp=mid.n();
   sf_warning("m2xp=%d n2xp=%d",m2xp, n2xp);

   fmidxp  = sf_floatalloc(m2xp*n2xp);
   ldataxp = sf_floatalloc(nxz*m2xp);
   rdataxp = sf_floatalloc(n2xp*nk);

   map2d1d(fmidxp, mid, m2xp, n2xp);

   iC ( samplexp(md,lid,mat) );
   map2d1d(ldataxp, mat, nxz, m2xp);

   iC ( samplexp(rid,nd,mat) );
   map2d1d(rdataxp, mat, n2xp, nk);

   /********* reconsturct W-matrix for checking **********/
   float **w;
   std::valarray<float> wx1(nk), wx2(nk);
   float **ap;
   if(ireconstruct==1)
   {
      w = sf_floatalloc2(nk,nxz);

      sf_warning("reconstruct x-component operators based on low-rank approx.");
      reconstruct(w, ldataxp, fmidxp, rdataxp, nxz, nk, m2xp, n2xp);

      /* error check between accurate and low-rank approx. */
      float sumall=0.0;
      ap=sf_floatalloc2(nz, nx);

      oRSF Errxp1("Errxp1");
      oRSF Errxp2("Errxp2");

      Errxp1.put("n1",nkz);
      Errxp1.put("n2",nkx);
      Errxp1.put("d1",dkz);
      Errxp1.put("d2",dkx);
      Errxp1.put("o1",kz0);
      Errxp1.put("o2",kx0);

      Errxp1.put("label1","kz");
      Errxp1.put("label2","kx");
      Errxp1.put("unit1","2*pi/m");
      Errxp1.put("unit2","2*pi/m");

      Errxp2.put("n1",nkz);
      Errxp2.put("n2",nkx);
      Errxp2.put("d1",dkz);
      Errxp2.put("d2",dkx);
      Errxp2.put("o1",kz0);
      Errxp2.put("o2",kx0);

      Errxp2.put("label1","kz");
      Errxp2.put("label2","kx");
      Errxp2.put("unit1","2*pi/m");
      Errxp2.put("unit2","2*pi/m");

      for(int im=0;im<nxz;im++){
        if(im%200==0) sf_warning("---------------------------------------- finish %d percent ---",(int)(im*100.0/nxz));

        //calculate P-wave's polarization operators directly for X-component
        polxp2dvti(ap, nx, nz, im);

        float sum=0.0;
        int k=0;
        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
        {
           float err=w[im][k]-ap[i][j];
           //if(fabs(err/ap[i][j])>0.01) sf_warning("Relative Error > 1%:  im=%d i=%d j=%d ikx=%d ikz=%d w=%f a=%f",im,i,j,k/nkz,k%nkz,w[im][k], ap[i][j]);
           sum += err*err;
           if(im==(nx/2-1)*nz+nz/4-1) wx1[k] = err;
           if(im==(nx/2-1)*nz+(nz*3)/4-1) wx2[k] = err;
           k++;
         }
         if(im%200==0) sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);
         sumall += sum/nxz;
       }// im loop
       sf_warning("X-component Low-rank average L2-norm error: %20.18f",sumall);

       Errxp1 << wx1;
       Errxp2 << wx2;

      oRSF Polxp1("Polxp1");
      oRSF Polxp2("Polxp2");

      Polxp1.put("n1",nkz);
      Polxp1.put("n2",nkx);
      Polxp1.put("d1",dkz);
      Polxp1.put("d2",dkx);
      Polxp1.put("o1",kz0);
      Polxp1.put("o2",kx0);

      Polxp1.put("label1","kz");
      Polxp1.put("label2","kx");
      Polxp1.put("unit1","2*pi/m");
      Polxp1.put("unit2","2*pi/m");

      Polxp2.put("n1",nkz);
      Polxp2.put("n2",nkx);
      Polxp2.put("d1",dkz);
      Polxp2.put("d2",dkx);
      Polxp2.put("o1",kz0);
      Polxp2.put("o2",kx0);

      Polxp2.put("label1","kz");
      Polxp2.put("label2","kx");
      Polxp2.put("unit1","2*pi/m");
      Polxp2.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++){
         wx1[i]=w[(nx/2-1)*nz+nz/4-1][i]; 
         wx2[i]=w[(nx/2-1)*nz+(nz*3)/4-1][i]; 
      }

      Polxp1 << wx1;
      Polxp2 << wx2;
   }//ireconstruct==1

   /********* low rank decomposition p-wave, z-component **********/
   iC( ddlowrank(nxz,nk,samplezp,eps,npk,lid,rid,mid) );
   m2zp=mid.m();
   n2zp=mid.n();
   sf_warning("m2zp=%d n2zp=%d",m2zp, n2zp);

   fmidzp  = sf_floatalloc(m2zp*n2zp);
   ldatazp = sf_floatalloc(nxz*m2zp);
   rdatazp = sf_floatalloc(n2zp*nk);

   map2d1d(fmidzp, mid, m2zp, n2zp);

   iC ( samplezp(md,lid,mat) );
   map2d1d(ldatazp, mat, nxz, m2zp);

   iC ( samplezp(rid,nd,mat) );
   map2d1d(rdatazp, mat, n2zp, nk);

   /********* reconsturct W-matrix for checking **********/
   std::valarray<float> wz1(nk), wz2(nk);
   if(ireconstruct==1)
   {
      w = sf_floatalloc2(nk,nxz);

      sf_warning("reconstruct z-component operators based on low-rank approx.");
      reconstruct(w, ldatazp, fmidzp, rdatazp, nxz, nk, m2zp, n2zp);

      /* error check between accurate and low-rank approx. */
      float sumall=0.0;
      ap=sf_floatalloc2(nz, nx);

      oRSF Errzp1("Errzp1");
      oRSF Errzp2("Errzp2");

      Errzp1.put("n1",nkz);
      Errzp1.put("n2",nkx);
      Errzp1.put("d1",dkz);
      Errzp1.put("d2",dkx);
      Errzp1.put("o1",kz0);
      Errzp1.put("o2",kx0);

      Errzp1.put("label1","kz");
      Errzp1.put("label2","kx");
      Errzp1.put("unit1","2*pi/m");
      Errzp1.put("unit2","2*pi/m");

      Errzp2.put("n1",nkz);
      Errzp2.put("n2",nkx);
      Errzp2.put("d1",dkz);
      Errzp2.put("d2",dkx);
      Errzp2.put("o1",kz0);
      Errzp2.put("o2",kx0);

      Errzp2.put("label1","kz");
      Errzp2.put("label2","kx");
      Errzp2.put("unit1","2*pi/m");
      Errzp2.put("unit2","2*pi/m");

      for(int im=0;im<nxz;im++){
        if(im%200==0) sf_warning("---------------------------------------- finish %d percent ---",(int)(im*100.0/nxz));

        //calculate P-wave's polarization operators directly for X-component
        polzp2dvti(ap, nx, nz, im);

        float sum=0.0;
        int k=0;
        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
        {
           float err=w[im][k]-ap[i][j];
           //if(fabs(err/ap[i][j])>0.01) sf_warning("Relative Error > 1%:  im=%d i=%d j=%d ikx=%d ikz=%d w=%f a=%f",im,i,j,k/nkz,k%nkz,w[im][k], ap[i][j]);
           sum += err*err;
           if(im==(nx/2-1)*nz+nz/4-1) wz1[k] = err;
           if(im==(nx/2-1)*nz+(nz*3)/4-1) wz2[k] = err;
           k++;
         }
         if(im%200==0) sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);
         sumall += sum/nxz;
       }
       // im loop
       sf_warning("Z-component Low-rank average L2-norm error: %20.18f",sumall);

      Errzp1 << wz1;
      Errzp2 << wz2;

      oRSF Polzp1("Polzp1");
      oRSF Polzp2("Polzp2");

      Polzp1.put("n1",nkz);
      Polzp1.put("n2",nkx);
      Polzp1.put("d1",dkz);
      Polzp1.put("d2",dkx);
      Polzp1.put("o1",kz0);
      Polzp1.put("o2",kx0);

      Polzp1.put("label1","kz");
      Polzp1.put("label2","kx");
      Polzp1.put("unit1","2*pi/m");
      Polzp1.put("unit2","2*pi/m");

      Polzp2.put("n1",nkz);
      Polzp2.put("n2",nkx);
      Polzp2.put("d1",dkz);
      Polzp2.put("d2",dkx);
      Polzp2.put("o1",kz0);
      Polzp2.put("o2",kx0);

      Polzp2.put("label1","kz");
      Polzp2.put("label2","kx");
      Polzp2.put("unit1","2*pi/m");
      Polzp2.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++){
         wz1[i]=w[(nx/2-1)*nz+nz/4-1][i]; 
         wz2[i]=w[(nx/2-1)*nz+(nz*3)/4-1][i]; 
      }
      
      Polzp1 << wz1;
      Polzp2 << wz2;
   }//ireconstruct==1

   /*****************************************************************************
   *  Calculating SV-wave polarization deviation operator for wave-mode separation
   * ***************************************************************************/
   /********* low rank decomposition SV-wave, x-component **********/
   int   m2xs, n2xs, m2zs, n2zs;
   float *ldataxs, *fmidxs, *rdataxs;
   float *ldatazs, *fmidzs, *rdatazs;

   iC( ddlowrank(nxz,nk,samplexs,eps,npk,lid,rid,mid) );
   m2xs=mid.m();
   n2xs=mid.n();
   sf_warning("m2xs=%d n2xs=%d",m2xs, n2xs);

   fmidxs  = sf_floatalloc(m2xs*n2xs);
   ldataxs = sf_floatalloc(nxz*m2xs);
   rdataxs = sf_floatalloc(n2xs*nk);

   map2d1d(fmidxs, mid, m2xs, n2xs);

   iC ( samplexs(md,lid,mat) );
   map2d1d(ldataxs, mat, nxz, m2xs);

   iC ( samplexs(rid,nd,mat) );
   map2d1d(rdataxs, mat, n2xs, nk);

   /********* reconsturct W-matrix for checking **********/
   if(ireconstruct==1)
   {
      w = sf_floatalloc2(nk,nxz);

      sf_warning("reconstruct x-component operators based on low-rank approx.");
      reconstruct(w, ldataxs, fmidxs, rdataxs, nxz, nk, m2xs, n2xs);

      /* error check between accurate and low-rank approx. */
      float sumall=0.0;
      ap=sf_floatalloc2(nz, nx);

      oRSF Errxs1("Errxs1");
      oRSF Errxs2("Errxs2");

      Errxs1.put("n1",nkz);
      Errxs1.put("n2",nkx);
      Errxs1.put("d1",dkz);
      Errxs1.put("d2",dkx);
      Errxs1.put("o1",kz0);
      Errxs1.put("o2",kx0);

      Errxs1.put("label1","kz");
      Errxs1.put("label2","kx");
      Errxs1.put("unit1","2*pi/m");
      Errxs1.put("unit2","2*pi/m");

      Errxs2.put("n1",nkz);
      Errxs2.put("n2",nkx);
      Errxs2.put("d1",dkz);
      Errxs2.put("d2",dkx);
      Errxs2.put("o1",kz0);
      Errxs2.put("o2",kx0);

      Errxs2.put("label1","kz");
      Errxs2.put("label2","kx");
      Errxs2.put("unit1","2*pi/m");
      Errxs2.put("unit2","2*pi/m");

      for(int im=0;im<nxz;im++){
        if(im%200==0) sf_warning("---------------------------------------- finish %d percent ---",(int)(im*100.0/nxz));

        //calculate P-wave's polarization operators directly for X-component
        polxs2dvti(ap, nx, nz, im);

        float sum=0.0;
        int k=0;
        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
        {
           float err=w[im][k]-ap[i][j];
           //if(fabs(err/ap[i][j])>0.01) sf_warning("Relative Error > 1%:  im=%d i=%d j=%d ikx=%d ikz=%d w=%f a=%f",im,i,j,k/nkz,k%nkz,w[im][k], ap[i][j]);
           if(im==(nx/2-1)*nz+nz/4-1) wx1[k] = err;
           if(im==(nx/2-1)*nz+(nz*3)/4-1) wx2[k] = err;
           sum += err*err;
           k++;
         }
         if(im%200==0) sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);
         sumall += sum/nxz;
       }
       // im loop
       sf_warning("X-component Low-rank average L2-norm error: %20.18f",sumall);

      Errxs1 << wx1;
      Errxs2 << wx2;

      oRSF Polxs1("Polxs1");
      oRSF Polxs2("Polxs2");

      Polxs1.put("n1",nkz);
      Polxs1.put("n2",nkx);
      Polxs1.put("d1",dkz);
      Polxs1.put("d2",dkx);
      Polxs1.put("o1",kz0);
      Polxs1.put("o2",kx0);

      Polxs1.put("label1","kz");
      Polxs1.put("label2","kx");
      Polxs1.put("unit1","2*pi/m");
      Polxs1.put("unit2","2*pi/m");

      Polxs2.put("n1",nkz);
      Polxs2.put("n2",nkx);
      Polxs2.put("d1",dkz);
      Polxs2.put("d2",dkx);
      Polxs2.put("o1",kz0);
      Polxs2.put("o2",kx0);

      Polxs2.put("label1","kz");
      Polxs2.put("label2","kx");
      Polxs2.put("unit1","2*pi/m");
      Polxs2.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++){
         wx1[i]=w[(nx/2-1)*nz+nz/4-1][i]; 
         wx2[i]=w[(nx/2-1)*nz+(nz*3)/4-1][i]; 
      }
      Polxs1 << wx1;
      Polxs2 << wx2;
      
   }//ireconstruct==1

   /********* low rank decomposition SV-wave, z-component **********/
   iC( ddlowrank(nxz,nk,samplezs,eps,npk,lid,rid,mid) );
   m2zs=mid.m();
   n2zs=mid.n();
   sf_warning("m2zs=%d n2zs=%d",m2zs, n2zs);

   fmidzs  = sf_floatalloc(m2zs*n2zs);
   ldatazs = sf_floatalloc(nxz*m2zs);
   rdatazs = sf_floatalloc(n2zs*nk);

   map2d1d(fmidzs, mid, m2zs, n2zs);

   iC ( samplezs(md,lid,mat) );
   map2d1d(ldatazs, mat, nxz, m2zs);

   iC ( samplezs(rid,nd,mat) );
   map2d1d(rdatazs, mat, n2zs, nk);

   /********* reconsturct W-matrix for checking **********/
   if(ireconstruct==1)
   {
      w = sf_floatalloc2(nk,nxz);

      sf_warning("reconstruct z-component operators based on low-rank approx.");
      reconstruct(w, ldatazs, fmidzs, rdatazs, nxz, nk, m2zs, n2zs);

      /* error check between accurate and low-rank approx. */
      float sumall=0.0;
      ap=sf_floatalloc2(nz, nx);

      oRSF Errzs1("Errzs1");
      oRSF Errzs2("Errzs2");

      Errzs1.put("n1",nkz);
      Errzs1.put("n2",nkx);
      Errzs1.put("d1",dkz);
      Errzs1.put("d2",dkx);
      Errzs1.put("o1",kz0);
      Errzs1.put("o2",kx0);

      Errzs1.put("label1","kz");
      Errzs1.put("label2","kx");
      Errzs1.put("unit1","2*pi/m");
      Errzs1.put("unit2","2*pi/m");

      Errzs2.put("n1",nkz);
      Errzs2.put("n2",nkx);
      Errzs2.put("d1",dkz);
      Errzs2.put("d2",dkx);
      Errzs2.put("o1",kz0);
      Errzs2.put("o2",kx0);

      Errzs2.put("label1","kz");
      Errzs2.put("label2","kx");
      Errzs2.put("unit1","2*pi/m");
      Errzs2.put("unit2","2*pi/m");

      for(int im=0;im<nxz;im++){
        if(im%200==0) sf_warning("---------------------------------------- finish %d percent ---",(int)(im*100.0/nxz));

        //calculate P-wave's polarization operators directly for X-component
        polzs2dvti(ap, nx, nz, im);

        float sum=0.0;
        int k=0;
        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
        {
           float err=w[im][k]-ap[i][j];
           //if(fabs(err/ap[i][j])>0.01) sf_warning("Relative Error > 1%:  im=%d i=%d j=%d ikx=%d ikz=%d w=%f a=%f",im,i,j,k/nkz,k%nkz,w[im][k], ap[i][j]);
           sum += err*err;
           if(im==(nx/2-1)*nz+nz/4-1) wz1[k] = err;
           if(im==(nx/2-1)*nz+(nz*3)/4-1) wz2[k] = err;
           k++;
         }
         if(im%200==0) sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);
         sumall += sum/nxz;
       }
       // im loop
       sf_warning("Z-component Low-rank average L2-norm error: %20.18f",sumall);

      Errzs1 << wz1;
      Errzs2 << wz2;

      oRSF Polzs1("Polzs1");
      oRSF Polzs2("Polzs2");

      Polzs1.put("n1",nkz);
      Polzs1.put("n2",nkx);
      Polzs1.put("d1",dkz);
      Polzs1.put("d2",dkx);
      Polzs1.put("o1",kz0);
      Polzs1.put("o2",kx0);

      Polzs1.put("label1","kz");
      Polzs1.put("label2","kx");
      Polzs1.put("unit1","2*pi/m");
      Polzs1.put("unit2","2*pi/m");

      Polzs2.put("n1",nkz);
      Polzs2.put("n2",nkx);
      Polzs2.put("d1",dkz);
      Polzs2.put("d2",dkx);
      Polzs2.put("o1",kz0);
      Polzs2.put("o2",kx0);

      Polzs2.put("label1","kz");
      Polzs2.put("label2","kx");
      Polzs2.put("unit1","2*pi/m");
      Polzs2.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++){
         wz1[i]=w[(nx/2-1)*nz+nz/4-1][i]; 
         wz2[i]=w[(nx/2-1)*nz+(nz*3)/4-1][i]; 
      }
      Polzs1 << wz1;
      Polzs2 << wz2;

      free(*ap);
      free(*w);
   }//ireconstruct==1

   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************End of Calculating Projection Deviation Operator****************/

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   /*  wavelet parameter for source definition */
   float A, f0, t0;
   f0=30.0;                  
   t0=0.04;                  
   A=1;                  

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

   /* setup I/O files */
   oRSF Elasticx("out"),Elasticz("Elasticz");
   oRSF ElasticSepP("ElasticSepP");
   oRSF ElasticSepSV("ElasticSepSV");

   Elasticx.put("n1",nkz);
   Elasticx.put("n2",nkx);
   Elasticx.put("d1",dz/1000);
   Elasticx.put("d2",dx/1000);
   Elasticx.put("o1",fz/1000);
   Elasticx.put("o2",fx/1000);

   Elasticz.put("n1",nkz);
   Elasticz.put("n2",nkx);
   Elasticz.put("d1",dz/1000);
   Elasticz.put("d2",dx/1000);
   Elasticz.put("o1",fz/1000);
   Elasticz.put("o2",fx/1000);

   ElasticSepP.put("n1",nkz);
   ElasticSepP.put("n2",nkx);
   ElasticSepP.put("d1",dz/1000);
   ElasticSepP.put("d2",dx/1000);
   ElasticSepP.put("o1",fz/1000);
   ElasticSepP.put("o2",fx/1000);

   ElasticSepSV.put("n1",nkz);
   ElasticSepSV.put("n2",nkx);
   ElasticSepSV.put("d1",dz/1000);
   ElasticSepSV.put("d2",dx/1000);
   ElasticSepSV.put("o1",fz/1000);
   ElasticSepSV.put("o2",fx/1000);

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

    std::valarray<float> x(nxz);
    std::valarray<float> y(nxz);
    std::valarray<float> z(nxz);

    float **vpp, **vss, **epp, **dee;

    vpp=sf_floatalloc2(nz,nx);
    vss=sf_floatalloc2(nz,nx);
    epp=sf_floatalloc2(nz,nx);
    dee=sf_floatalloc2(nz,nx);

    k=0;
    for(i=0;i<nx;i++)
    for(j=0;j<nz;j++)
    {
       vpp[i][j]=vp[k];
       vss[i][j]=vs[k];
       epp[i][j]=ep[k];
       dee[i][j]=de[k];
       k++;
    }

    t4=clock();
    timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
    sf_warning("CPU time for preparing for wavefield modeling: %f(second)",timespent);

    float *pp, *qq;
    pp=sf_floatalloc(nxz);
    qq=sf_floatalloc(nxz);

    int ii, jj, im, jm;

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

    int iflag=0;

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
        /* fwpvtielastic: forward-propagating using original elastic equation of displacement in VTI media*/
        fwpvtielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
                      dx, dz, nx, nz, nxpad, nzpad, vpp, vss, epp, dee);

        /******* output wavefields: component and divergence *******/
        if(it==ns-1)
	{
              k=0;
	      for(i=0;i<nx;i++)
              {
                   im=i+M;
		   for(j=0;j<nz;j++)
		   {
                       jm=j+M;

                       x[k] = pp[k] = p3[im][jm];
                       y[k] = qq[k] = q3[im][jm];

                       k++;      
		    }
              }// i loop
              Elasticx<<x;
              Elasticz<<y;

              t44=clock();
              timespent=(float)(t44-t4)/CLOCKS_PER_SEC;
              sf_warning("CPU time for wavefield modeling: %f(second)",timespent);

              /* separate qP wave  */
              sf_warning("separate qP-wave based on lowrank decomp."); 
              seplowrank2d(ldataxp,rdataxp,fmidxp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xp,n2xp,iflag);
              seplowrank2d(ldatazp,rdatazp,fmidzp,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zp,n2zp,iflag);

              for(i=0;i<nxz;i++)
		  z[i]=pp[i]+qq[i];

              ElasticSepP<<z;
          
              k=0;
	      for(i=0;i<nx;i++)
              {
                   im=i+M;
		   for(j=0;j<nz;j++)
		   {
                       jm=j+M;

                       pp[k] = p3[im][jm];
                       qq[k] = q3[im][jm];

                       k++;      
		    }
              }// i loop

              /* separate qSV wave  */
              sf_warning("separate qSV-wave based on lowrank decomp."); 
              seplowrank2d(ldataxs,rdataxs,fmidxs,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xs,n2xs,iflag);
              seplowrank2d(ldatazs,rdatazs,fmidzs,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zs,n2zs,iflag);

              for(i=0;i<nxz;i++)
		  z[i]=pp[i]+qq[i];

              ElasticSepSV<<z;
          
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
    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

    free(ldataxp);
    free(ldatazp);
    free(rdataxp);
    free(rdatazp);
    free(fmidxp);
    free(fmidzp);
    free(ldataxs);
    free(ldatazs);
    free(rdataxs);
    free(rdatazs);
    free(fmidxs);
    free(fmidzs);

    free(pp);
    free(qq);

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

    free(coeff_2dx);
    free(coeff_2dz);
    free(coeff_1dx);
    free(coeff_1dz);

    free(ijkx);
    free(ijkz);
exit(0);
}
/* P-wave dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            a11= c11*s*s+c44*c*c;
            a12= c13c44*s*c;
            a22= c44*s*s+c33*c*c;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*s + u2*c <0) {
               u1 = -u1;
               u2 = -u2;
            }

            resx(a,b) = u1;
              
         }// b loop
    }// a loop

    return 0;
}

/* P-wave dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezp(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];

            if(s==0&&c==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            a11= c11*s*s+c44*c*c;
            a12= c13c44*s*c;
            a22= c44*s*s+c33*c*c;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*s + u2*c <0) {
               u1 = -u1;
               u2 = -u2;
            }

            resz(a,b) = u2;
              
         }// b loop
    }// a loop

    return 0;
}

/* SV-wave dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexs(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            a11= c11*s*s+c44*c*c;
            a12= c13c44*s*c;
            a22= c44*s*s+c33*c*c;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[1][0];
            double u2=ve[1][1];

            /* get the closest direction to k */
            if(u1*c - u2*s <0) {
               u1 = -u1;
               u2 = -u2;
            }

            resx(a,b) = u1;
              
         }// b loop
    }// a loop

    return 0;
}

/* SV-wave dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezs(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

    double   aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/

    double c33, c44, c11, c13c44, a11, a12, a22;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];

            if(s==0&&c==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            a11= c11*s*s+c44*c*c;
            a12= c13c44*s*c;
            a22= c44*s*s+c33*c*c;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            double u1=ve[1][0];
            double u2=ve[1][1];

            /* get the closest direction to k */
            if(u1*c - u2*s <0) {
               u1 = -u1;
               u2 = -u2;
            }

            resz(a,b) = u2;
              
         }// b loop
    }// a loop

    return 0;
}

static void polxp2dvti(float **ap, int nx, int nz, int im)
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

                a11= c11*sx*sx + c44*cx*cx;
                a12= c13c44*sx*cx;
                a22= c44*sx*sx + c33*cx*cx;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[0][0];
               u2=ve[0][1];
               /* get the closest direction to k */
               if(u1*sx + u2*cx <0) {
                  u1 = -u1;
               }
	       ap[i][j]=(float)(u1);
          } /* j loop */
}

static void polzp2dvti(float **ap, int nx, int nz, int im)
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

                a11= c11*sx*sx + c44*cx*cx;
                a12= c13c44*sx*cx;
                a22= c44*sx*sx + c33*cx*cx;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[0][0];
               u2=ve[0][1];
               /* get the closest direction to k */
               if(u1*sx + u2*cx <0) {
                  u2 = -u2;
               }
	       ap[i][j]=(float)(u2);
          } /* j loop */
}

static void polxs2dvti(float **ap, int nx, int nz, int im)
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

                a11= c11*sx*sx + c44*cx*cx;
                a12= c13c44*sx*cx;
                a22= c44*sx*sx + c33*cx*cx;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[1][0];
               u2=ve[1][1];
               /* get the closest direction to k */
               if(u1*cx - u2*sx <0) {
                  u1 = -u1;
               }
	       ap[i][j]=(float)(u1);
          } /* j loop */
}

static void polzs2dvti(float **ap, int nx, int nz, int im)
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

                a11= c11*sx*sx + c44*cx*cx;
                a12= c13c44*sx*cx;
                a22= c44*sx*sx + c33*cx*cx;

                a[0][0] = a11;
                a[0][1] = a12;
                a[1][0] = a12;
                a[1][1] = a22;

               dsolveSymmetric22(a, ve, va);

               u1=ve[1][0];
               u2=ve[1][1];
               /* get the closest direction to k */
               if(u1*cx - u2*sx <0) {
                  u2 = -u2;
               }
	       ap[i][j]=(float)(u2);
          } /* j loop */
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

