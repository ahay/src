/* 3-D three-components wavefield modeling based on pseudo-pure P-wave equation
   and P-S separation using low-rank symbol approximation in 3D ORT media.

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
#include "seplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, ga, th, al;

static std::valarray<double> rkx, rky, rkz, rkk;

/* for low-rank decomp. */
static int samplexp3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int sampleyp3(vector<int>& rs, vector<int>& cs, DblNumMat& resy);
static int samplezp3(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static int samplexsh3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int sampleysh3(vector<int>& rs, vector<int>& cs, DblNumMat& resy);
static int samplezsh3(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static int samplexsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int sampleysv3(vector<int>& rs, vector<int>& cs, DblNumMat& resy);
static int samplezsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

/* conventional computation of projection deviation operator */
static void dev3dort(float ***ap, int nx, int ny, int nz, int im);

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
   par.get("ireconstruct",ireconstruct);
   par.get("iflagvti",iflagvti);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("ireconstruct=%d ",ireconstruct);
   sf_warning("iflagvti=%d ",iflagvti);
   sf_warning("read velocity model parameters");

   /* setup I/O files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), gam("gam"), the("the"), alph("alph");
   iRSF Fx("Elasticx"), Fy("Elasticy"),Fz("Elasticz");
   oRSF Fp("ElasticSepP"), Fsv("ElasticSepSV"),Fsh("ElasticSepSH");

   float az, ax, ay;

   /* Read/Write axes from Model*/
   int nxv, nyv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("n3",nyv);
   vp0.get("o1",az);
   vp0.get("o2",ax);
   vp0.get("o3",ay);
   float fxv, fyv, fzv;
   fxv=ax*1000.0;
   fyv=ay*1000.0;
   fzv=az*1000.0;
   float dxv, dyv, dzv;
   vp0.get("d1",az);
   vp0.get("d2",ax);
   vp0.get("d3",ay);
   dzv = az*1000.0;
   dxv = ax*1000.0;
   dyv = ay*1000.0;

   /* Read/Write axes from Wavefields*/
   int nx, ny, nz;
   Fx.get("n1",nz);
   Fx.get("n2",nx);
   Fx.get("n3",ny);
   if(nx!=nxv || ny!=nyv || nz!=nzv){
     sf_warning("Dimension not match between model and data !");
     sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
     sf_warning("nxv=%d nyv=%d nzv=%d ",nxv,nyv,nzv);
     exit(0);
   }
   Fx.get("o1",az);
   Fx.get("o2",ax);
   Fx.get("o3",ay);
   float fx, fy, fz;
   fx=ax*1000.0;
   fy=ay*1000.0;
   fz=az*1000.0;
   if(fabs(fx-fxv)>0.1 || fabs(fy-fyv)>0.1 || fabs(fz-fzv)>0.1){
     sf_warning("Coorinate original point not match between model and data !");
     sf_warning("fx=%d fy=%d fz=%d ",fx,fy,fz);
     sf_warning("fxv=%d fyv=%d fzv=%d ",fxv,fyv,fzv);
   }
   float dx, dy, dz;
   Fx.get("d1",az);
   Fx.get("d2",ax);
   Fx.get("d3",ay);
   dz = az*1000.0;
   dx = ax*1000.0;
   dy = ay*1000.0;
   if(fabs(dx-dxv)>0.1 || fabs(dy-dyv)>0.1 || fabs(dz-dzv)>0.1){
     sf_warning("Sampling step not match between model and data !");
     sf_warning("dx=%d dy=%d dz=%d ",dx,dy,dz);
     sf_warning("dxv=%d dyv=%d dzv=%d ",dxv,dyv,dzv);
   }

   sf_warning("fx=%f fy=%f fz=%f dx=%f dy=%f dz=%f",fx,fy,fz,dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d ", nx,ny,nz);

   /* wave modeling space */
   int nxyz=nx*ny*nz;

   vp.resize(nxyz);
   vs.resize(nxyz);
   ep.resize(nxyz);
   de.resize(nxyz);
   ga.resize(nxyz);
   th.resize(nxyz);
   al.resize(nxyz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;
   gam>>ga;
   the>>th;
   alph>>al;
  
   for(int i=0;i<nxyz;i++){
      th[i] *= PI/180.0;
      al[i] *= PI/180.0;
   } 

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nxyz=%d nk=%d",nxyz,nk);

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
                    low rank decomposition for P-wave's operators
   * ***************************************************************************/
   int   m2yp, n2yp, m2xp, n2xp, m2zp, n2zp;

   /********* low rank decomposition p-wave, y-component **********/
   iC( ddlowrank(nxyz,nk,sampleyp3,eps,npk,lid,rid,mid) );
   m2yp=mid.m();
   n2yp=mid.n();
   sf_warning("lowrank-p-y:m2yp=%d n2yp=%d",m2yp, n2yp);

   float *ldataxp, *fmidxp, *rdataxp;
   float *ldatayp, *fmidyp, *rdatayp;
   float *ldatazp, *fmidzp, *rdatazp;

   fmidyp  = sf_floatalloc(m2yp*n2yp);
   ldatayp = sf_floatalloc(nxyz*m2yp);
   rdatayp = sf_floatalloc(n2yp*nk);

   map2d1d(fmidyp, mid, m2yp, n2yp);

   iC ( sampleyp3(md,lid,mat) );
   map2d1d(ldatayp, mat, nxyz, m2yp);

   iC ( sampleyp3(rid,nd,mat) );
   map2d1d(rdatayp, mat, n2yp, nk);

   /********* reconsturct W-matrix for checking **********/
   float **w;
   float ***ap;
   std::valarray<float> wy(nk);
   if(ireconstruct==1)
   {
      w = sf_floatalloc2(nk,nxyz);

      sf_warning("reconstruct y-component operators based on low-rank approx.");
      reconstruct(w, ldatayp, fmidyp, rdatayp, nxyz, nk, m2yp, n2yp);

      /* error check between accurate and low-rank approx. */
      ap=sf_floatalloc3(nz, nx, ny);

      float sumall=0.0;
      //for(int im=0;im<nxyz;im++){
      for(int im=nxyz/2;im<=nxyz/2;im++){
        if(im%200==0||im==nxyz/2) sf_warning("---------------------------------------- finish %d percent ---",(int)(im*100.0/nxyz));

        //calculate P-wave's polarization operators directly for Y-component
        dev3dort(ap, nkx, nky, nkz, im);

        float sum=0.0;
        int l, k=0;
        for(l=0;l<nky;l++)
        for(i=0;i<nkx;i++)
        for(j=0;j<nkz;j++)
        {
           if((im%200==0||im==nxyz/2)&&l==0&&i==0&&j%50==0) sf_warning("w=%20.16f ap=%20.16f",w[im][k],ap[l][i][j]);
           float err=w[im][k]-ap[l][i][j];
           sum += err*err;
           if(im==nxyz/2) wy[k] = err*1000000;
           k++;
         }
         if(im%200==0||im==nxyz/2) sf_warning("Y-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxyz);
         sumall += sum/nxyz;
       }// im loop
       sf_warning("Y-component Low-rank average L2-norm error: %20.18f",sumall);

      oRSF Polyp("out"), Erryp("Erryp");

      Erryp.put("n1",nkz);
      Erryp.put("n2",nkx);
      Erryp.put("n3",nky);
      Erryp.put("d1",dkz);
      Erryp.put("d2",dkx);
      Erryp.put("d3",dky);
      Erryp.put("o1",kz0);
      Erryp.put("o2",kx0);
      Erryp.put("o3",ky0);

      Erryp.put("label1","kz");
      Erryp.put("label2","kx");
      Erryp.put("label3","ky");
      Erryp.put("unit1","2*pi/m");
      Erryp.put("unit2","2*pi/m");
      Erryp.put("unit3","2*pi/m");

      Erryp << wy;

      Polyp.put("n1",nkz);
      Polyp.put("n2",nkx);
      Polyp.put("n3",nky);
      Polyp.put("d1",dkz);
      Polyp.put("d2",dkx);
      Polyp.put("d3",dky);
      Polyp.put("o1",kz0);
      Polyp.put("o2",kx0);
      Polyp.put("o3",ky0);

      Polyp.put("label1","kz");
      Polyp.put("label2","kx");
      Polyp.put("label3","ky");
      Polyp.put("unit1","2*pi/m");
      Polyp.put("unit2","2*pi/m");
      Polyp.put("unit3","2*pi/m");

      for(i=0;i<nk;i++)
         wy[i]=w[nxyz/2][i];

      Polyp << wy;

      free(**ap);

   }//ireconstruct==1
   /********* low rank decomposition p-wave, x-component **********/
   iC( ddlowrank(nxyz,nk,samplexp3,eps,npk,lid,rid,mid) );
   m2xp=mid.m();
   n2xp=mid.n();
   sf_warning("lowrank-p-x:m2xp=%d n2xp=%d",m2xp, n2xp);

   fmidxp  = sf_floatalloc(m2xp*n2xp);
   ldataxp = sf_floatalloc(nxyz*m2xp);
   rdataxp = sf_floatalloc(n2xp*nk);

   map2d1d(fmidxp, mid, m2xp, n2xp);

   iC ( samplexp3(md,lid,mat) );
   map2d1d(ldataxp, mat, nxyz, m2xp);

   iC ( samplexp3(rid,nd,mat) );
   map2d1d(rdataxp, mat, n2xp, nk);

   /********* low rank decomposition p-wave, z-component **********/
   iC( ddlowrank(nxyz,nk,samplezp3,eps,npk,lid,rid,mid) );
   m2zp=mid.m();
   n2zp=mid.n();
   sf_warning("lowrank-p-z:m2zp=%d n2zp=%d",m2zp, n2zp);

   fmidzp  = sf_floatalloc(m2zp*n2zp);
   ldatazp = sf_floatalloc(nxyz*m2zp);
   rdatazp = sf_floatalloc(n2zp*nk);

   map2d1d(fmidzp, mid, m2zp, n2zp);

   iC ( samplezp3(md,lid,mat) );
   map2d1d(ldatazp, mat, nxyz, m2zp);

   iC ( samplezp3(rid,nd,mat) );
   map2d1d(rdatazp, mat, n2zp, nk);

   /*****************************************************************************
                    low rank decomposition for SV-wave's operators
   * ***************************************************************************/
   int   m2ysv, n2ysv, m2xsv, n2xsv, m2zsv, n2zsv;
   float *ldataxsv, *fmidxsv, *rdataxsv;
   float *ldataysv, *fmidysv, *rdataysv;
   float *ldatazsv, *fmidzsv, *rdatazsv;

   /********* low rank decomposition SV-wave, y-component **********/
   iC( ddlowrank(nxyz,nk,sampleysv3,eps,npk,lid,rid,mid) );
   m2ysv=mid.m();
   n2ysv=mid.n();
   sf_warning("lowrank-sv-y:m2ysv=%d n2ysv=%d",m2ysv, n2ysv);

   fmidysv  = sf_floatalloc(m2ysv*n2ysv);
   ldataysv = sf_floatalloc(nxyz*m2ysv);
   rdataysv = sf_floatalloc(n2ysv*nk);

   map2d1d(fmidysv, mid, m2ysv, n2ysv);

   iC ( sampleysv3(md,lid,mat) );
   map2d1d(ldataysv, mat, nxyz, m2ysv);

   iC ( sampleysv3(rid,nd,mat) );
   map2d1d(rdataysv, mat, n2ysv, nk);

   /********* low rank decomposition SV-wave, x-component **********/

   iC( ddlowrank(nxyz,nk,samplexsv3,eps,npk,lid,rid,mid) );
   m2xsv=mid.m();
   n2xsv=mid.n();
   sf_warning("lowrank-sv-x:m2xsv=%d n2xsv=%d",m2xsv, n2xsv);

   fmidxsv  = sf_floatalloc(m2xsv*n2xsv);
   ldataxsv = sf_floatalloc(nxyz*m2xsv);
   rdataxsv = sf_floatalloc(n2xsv*nk);

   map2d1d(fmidxsv, mid, m2xsv, n2xsv);

   iC ( samplexsv3(md,lid,mat) );
   map2d1d(ldataxsv, mat, nxyz, m2xsv);

   iC ( samplexsv3(rid,nd,mat) );
   map2d1d(rdataxsv, mat, n2xsv, nk);

   /********* low rank decomposition SV-wave, z-component **********/

   iC( ddlowrank(nxyz,nk,samplezsv3,eps,npk,lid,rid,mid) );
   m2zsv=mid.m();
   n2zsv=mid.n();
   sf_warning("lowrank-sv-z:m2zsv=%d n2zsv=%d",m2zsv, n2zsv);

   fmidzsv  = sf_floatalloc(m2zsv*n2zsv);
   ldatazsv = sf_floatalloc(nxyz*m2zsv);
   rdatazsv = sf_floatalloc(n2zsv*nk);

   map2d1d(fmidzsv, mid, m2zsv, n2zsv);

   iC ( samplezsv3(md,lid,mat) );
   map2d1d(ldatazsv, mat, nxyz, m2zsv);

   iC ( samplezsv3(rid,nd,mat) );
   map2d1d(rdatazsv, mat, n2zsv, nk);
   
   /*****************************************************************************
                    low rank decomposition for SH-wave's operators
   * ***************************************************************************/
   int   m2ysh, n2ysh, m2xsh, n2xsh, m2zsh, n2zsh;
   float *ldataxsh, *fmidxsh, *rdataxsh;
   float *ldataysh, *fmidysh, *rdataysh;
   float *ldatazsh, *fmidzsh, *rdatazsh;

   /********* low rank decomposition SH-wave, y-component **********/
   iC( ddlowrank(nxyz,nk,sampleysh3,eps,npk,lid,rid,mid) );
   m2ysh=mid.m();
   n2ysh=mid.n();
   sf_warning("lowrank-sh-y:m2ysh=%d n2ysh=%d",m2ysh, n2ysh);

   fmidysh  = sf_floatalloc(m2ysh*n2ysh);
   ldataysh = sf_floatalloc(nxyz*m2ysh);
   rdataysh = sf_floatalloc(n2ysh*nk);

   map2d1d(fmidysh, mid, m2ysh, n2ysh);

   iC ( sampleysh3(md,lid,mat) );
   map2d1d(ldataysh, mat, nxyz, m2ysh);

   iC ( sampleysh3(rid,nd,mat) );
   map2d1d(rdataysh, mat, n2ysh, nk);

   /********* low rank decomposition SH-wave, x-component **********/
   iC( ddlowrank(nxyz,nk,samplexsh3,eps,npk,lid,rid,mid) );
   m2xsh=mid.m();
   n2xsh=mid.n();
   sf_warning("lowrank-sh-x:m2xsh=%d n2xsh=%d",m2xsh, n2xsh);

   fmidxsh  = sf_floatalloc(m2xsh*n2xsh);
   ldataxsh = sf_floatalloc(nxyz*m2xsh);
   rdataxsh = sf_floatalloc(n2xsh*nk);

   map2d1d(fmidxsh, mid, m2xsh, n2xsh);

   iC ( samplexsh3(md,lid,mat) );
   map2d1d(ldataxsh, mat, nxyz, m2xsh);

   iC ( samplexsh3(rid,nd,mat) );
   map2d1d(rdataxsh, mat, n2xsh, nk);

   /********* low rank decomposition SH-wave, z-component **********/
   if(iflagvti==0){
      iC( ddlowrank(nxyz,nk,samplezsh3,eps,npk,lid,rid,mid) );
      m2zsh=mid.m();
      n2zsh=mid.n();
      sf_warning("lowrank-sh-z:m2zsh=%d n2zsh=%d",m2zsh, n2zsh);

      fmidzsh  = sf_floatalloc(m2zsh*n2zsh);
      ldatazsh = sf_floatalloc(nxyz*m2zsh);
      rdatazsh = sf_floatalloc(n2zsh*nk);

      map2d1d(fmidzsh, mid, m2zsh, n2zsh);

      iC ( samplezsh3(md,lid,mat) );
      map2d1d(ldatazsh, mat, nxyz, m2zsh);

      iC ( samplezsh3(rid,nd,mat) );
      map2d1d(rdatazsh, mat, n2zsh, nk);
   }else{
      sf_warning("For 3D VTI Media, Z-comp. of SH's Polarization is always zero !!!!");
   }
   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);


   /****************begin to input wavefield (snapshot) and separation ****************/
   /****************begin to input wavefield (snapshot) and separation ****************/
   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

    float *px, *py, *pz;
    float *pp;
    std::valarray<float> p(nxyz);

    px=sf_floatalloc(nxyz);
    py=sf_floatalloc(nxyz);
    pz=sf_floatalloc(nxyz);
    pp=sf_floatalloc(nxyz);

    int iflag=0;

    Fx >> p;
    for(k=0;k<nxyz;k++) px[k] = p[k];
    Fy >> p;
    for(k=0;k<nxyz;k++) py[k] = p[k];
    Fz >> p;
    for(k=0;k<nxyz;k++) pz[k] = p[k];

        Fp.put("n1",nz);
        Fp.put("n2",nx);
        Fp.put("n3",ny);
        Fp.put("d1",dz);
        Fp.put("d2",dx);
        Fp.put("d3",dy);
        Fp.put("o1",fz);
        Fp.put("o2",fx);
        Fp.put("o3",fy);
        Fp.put("label1","z");
        Fp.put("label2","x");
        Fp.put("label3","y");
        Fp.put("unit1","km");
        Fp.put("unit2","km");
        Fp.put("unit3","km");

        Fsv.put("n1",nz);
        Fsv.put("n2",nx);
        Fsv.put("n3",ny);
        Fsv.put("d1",dz);
        Fsv.put("d2",dx);
        Fsv.put("d3",dy);
        Fsv.put("o1",fz);
        Fsv.put("o2",fx);
        Fsv.put("o3",fy);
        Fsv.put("label1","z");
        Fsv.put("label2","x");
        Fsv.put("label3","y");
        Fsv.put("unit1","km");
        Fsv.put("unit2","km");
        Fsv.put("unit3","km");

        Fsh.put("n1",nz);
        Fsh.put("n2",nx);
        Fsh.put("n3",ny);
        Fsh.put("d1",dz);
        Fsh.put("d2",dx);
        Fsh.put("d3",dy);
        Fsh.put("o1",fz);
        Fsh.put("o2",fx);
        Fsh.put("o3",fy);
        Fsh.put("label1","z");
        Fsh.put("label2","x");
        Fsh.put("label3","y");
        Fsh.put("unit1","km");
        Fsh.put("unit2","km");
        Fsh.put("unit3","km");

    //puthead3x(Fp, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    //puthead3x(Fsv, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    //puthead3x(Fsh, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

    // separate qP wave  
    sf_warning("separate qP-wave based on lowrank decomp."); 
    for(k=0;k<nxyz;k++) p[k] = 0.0;
    for(k=0;k<nxyz;k++) pp[k] = px[k];
    seplowrank3d(ldataxp,rdataxp,fmidxp,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2xp,n2xp,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = py[k];
    seplowrank3d(ldatayp,rdatayp,fmidyp,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2yp,n2yp,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = pz[k];
    seplowrank3d(ldatazp,rdatazp,fmidzp,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2zp,n2zp,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    Fp << p; 

    // separate qSV wave  
    sf_warning("separate qSV-wave based on lowrank decomp."); 
    for(k=0;k<nxyz;k++) p[k] = 0.0;
    for(k=0;k<nxyz;k++) pp[k] = px[k];
    seplowrank3d(ldataxsv,rdataxsv,fmidxsv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2xsv,n2xsv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = py[k];
    seplowrank3d(ldataysv,rdataysv,fmidysv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2ysv,n2ysv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = pz[k];
    seplowrank3d(ldatazsv,rdatazsv,fmidzsv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2zsv,n2zsv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    Fsv << p; 

    // separate qSH wave  
    sf_warning("separate qSH-wave based on lowrank decomp."); 
    for(k=0;k<nxyz;k++) p[k] = 0.0;
    for(k=0;k<nxyz;k++) pp[k] = px[k];
    seplowrank3d(ldataxsh,rdataxsh,fmidxsh,px,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2xsh,n2xsh,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = py[k];
    seplowrank3d(ldataysh,rdataysh,fmidysh,py,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2ysh,n2ysh,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    for(k=0;k<nxyz;k++) pp[k] = pz[k];
    seplowrank3d(ldatazsh,rdatazsh,fmidzsh,pz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2zsh,n2zsh,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

    Fsh << p; 

    free(px);
    free(py);
    free(pz);
    free(pp);

    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

    free(ldataxp);
    free(ldatayp);
    free(ldatazp);
    free(rdataxp);
    free(rdatayp);
    free(rdatazp);
    free(fmidxp);
    free(fmidyp);
    free(fmidzp);
    free(ldataxsv);
    free(ldataysv);
    free(ldatazsv);
    free(rdataxsv);
    free(rdataysv);
    free(rdatazsv);
    free(fmidxsv);
    free(fmidysv);
    free(fmidzsv);

    if(iflagvti==0){
      free(ldataxsh);
      free(ldataysh);
      free(ldatazsh);
      free(rdataxsh);
      free(rdataysh);
      free(rdatazsh);
      free(fmidxsh);
      free(fmidysh);
      free(fmidzsh);
    }

    sf_warning("-------sucessful ending --------");
    exit(0);
}

/* P-wave y-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int sampleyp3(vector<int>& rs, vector<int>& cs, DblNumMat& resy)
{
    int nr = rs.size();
    int nc = cs.size();

    resy.resize(nr,nc);

    setvalue(resy,0.0);
 
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

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resy(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            resy(a,b) = upy;
              
         }// b loop
    }// a loop

    return 0;
}

/* P-wave x-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

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

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            resx(a,b) = upx;
              
         }// b loop
    }// a loop

    return 0;
}

/* P-wave z-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezp3(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

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

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            resz(a,b) = upz;
              
         }// b loop
    }// a loop

    return 0;
}

/* SH-wave y-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int sampleysh3(vector<int>& rs, vector<int>& cs, DblNumMat& resy)
{
    int nr = rs.size();
    int nc = cs.size();

    resy.resize(nr,nc);

    setvalue(resy,0.0);

    double kx, ky, kz;
    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resy(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            usx /= rk;
            usy /= rk;
            usz /= rk;

            resy(a,b) = usy;
              
         }// b loop
    }// a loop

    return 0;
}

/* SH-wave x-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexsh3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double kx, ky, kz;
    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            usx /= rk;
            usy /= rk;
            usz /= rk;

            resx(a,b) = usx;
              
         }// b loop
    }// a loop

    return 0;
}

/* SH-wave z-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezsh3(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

    double kx, ky, kz;
    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
            kz = -sin_th*kx0                    + cos_th*kz0;

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            //sf_warning("ushx=%f ushy=%f ushz=%f",usx,usy,usz);

            usx /= rk;
            usy /= rk;
            usz /= rk;

            resz(a,b) = usz;
              
         }// b loop
    }// a loop

    return 0;
}

/* SV-wave y-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int sampleysv3(vector<int>& rs, vector<int>& cs, DblNumMat& resy)
{
    int nr = rs.size();
    int nc = cs.size();

    resy.resize(nr,nc);

    setvalue(resy,0.0);

    double kx, ky, kz;
    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resy(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            double usx= ky*xn*upy - kx*yn*upy + kz*xn*upz - kx*zn*upz;
            double usy= kz*yn*upz - ky*zn*upz + kx*yn*upx - ky*xn*upx;
            double usz= kx*zn*upx - kz*xn*upx + ky*zn*upy - kz*yn*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            //usx /= rk;
            usy /= rk;
            //usz /= rk;

            resy(a,b) = usy;
              
         }// b loop
    }// a loop

    return 0;
}

/* SV-wave x-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double kx, ky, kz;
    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            double usx= ky*xn*upy - kx*yn*upy + kz*xn*upz - kx*zn*upz;
            double usy= kz*yn*upz - ky*zn*upz + kx*yn*upx - ky*xn*upx;
            double usz= kx*zn*upx - kz*xn*upx + ky*zn*upy - kz*yn*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            usx /= rk;
            //usy /= rk;
            //usz /= rk;

            resx(a,b) = usx;
              
         }// b loop
    }// a loop

    return 0;
}
/* SV-wave z-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

    double kx, ky, kz;
    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        double cos_al=cos(al[i]);
        double sin_al=sin(al[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_al*sin_th;
        double yn= sin_al*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx0 = rkx[cs[b]];
            double ky0 = rky[cs[b]];
            double kz0 = rkz[cs[b]];
            if(kx0==0&&ky0==0&&kz0==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            /* tilted axis proccessing */
            kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
            ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

            double usx= ky*xn*upy - kx*yn*upy + kz*xn*upz - kx*zn*upz;
            double usy= kz*yn*upz - ky*zn*upz + kx*yn*upx - ky*xn*upx;
            double usz= kx*zn*upx - kz*xn*upx + ky*zn*upy - kz*yn*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            
            //usx /= rk;
            //usy /= rk;
            usz /= rk;

            resz(a,b) = usz;
              
         }// b loop
    }// a loop

    return 0;
}
/* P-wave y-component polarization operator */
static void dev3dort(float ***ap, int nx, int ny, int nz, int im)
{
    int     i, j, k, l;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    double kx, ky, kz;
    double upx, upy, upz;

    double vp2 = vp[im]*vp[im];
    double vs2 = vs[im]*vs[im];
    double ep2 = 1.0+2*ep[im];
    double de2 = 1.0+2*de[im];
    double ga2 = 1.0+2*ga[im];

    double cos_al=cos(al[im]);
    double sin_al=sin(al[im]);
    double cos_th=cos(th[im]);
    double sin_th=sin(th[im]);

    for( l=0; l<ny ; l++)
    for( i=0; i<nx ; i++)
    for( j=0; j<nz ; j++)
         ap[l][i][j]=0.0;

    k=0;
    for( l=0; l<ny ; l++)
    for( i=0; i<nx ; i++)
    for( j=0; j<nz ; j++)
    {
         double kx0 = rkx[k];
         double ky0 = rky[k];
         double kz0 = rkz[k];
         if(kx0==0&&ky0==0&&kz0==0)
         {
            ap[l][i][j]=0.0;
            continue;
         }
         k++;

         c33=vp2;
         c44=vs2;
         c11=ep2*c33;
         c66=ga2*c44;
         c13c44=sqrt((de2*c33-c44)*(c33-c44));
         c11c66=c11-c66;

         /* tilted axis proccessing */
         kx = cos_th*cos_al*kx0 - sin_al*ky0 + sin_th*cos_al*kz0;
         ky = cos_th*sin_al*kx0 + cos_al*ky0 + sin_th*sin_al*kz0;
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

         ap[l][i][j] = (float)upy;
      }

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

/*static void puthead3x(oRSF Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
{
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("n3",n3);
        Fo.put("d1",d1);
        Fo.put("d2",d2);
        Fo.put("d3",d3);
        Fo.put("o1",o1);
        Fo.put("o2",o2);
        Fo.put("o3",o3);
        Fo.put("label1","z");
        Fo.put("label2","x");
        Fo.put("label3","y");
        Fo.put("unit1","km");
        Fo.put("unit2","km");
        Fo.put("unit3","km");
}*/

