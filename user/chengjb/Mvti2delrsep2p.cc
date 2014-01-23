/* 2-D two-components wavefield modeling based on original elastic anisotropic displacement
  wave equation and P-SV separation using low-rank symbol approximation in 2D TTI media.

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
#include "puthead.h"
#include "kykxkztaper.h"
#include "eigen2x2.h"
#include "seplowrank.h"
}

static std::valarray<float> vp, vs, ep, de;

static std::valarray<double> sinx, cosx;

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplexs(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezs(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

static void polxp2dtti(float **ap, int nx, int nz, int im);
static void polzp2dtti(float **ap, int nx, int nz, int im);

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
   float dt=0.004, xrec1=1000.0, zrec1=1000.0, xrec2=2000.0, zrec2=2000.0;

   int ireconstruct;   // flag for reconstruct the W matrix or not
   par.get("ireconstruct",ireconstruct);
   par.get("xrec1",xrec1);
   par.get("zrec1",zrec1);
   par.get("xrec2",xrec2);
   par.get("zrec2",zrec2);

   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("ireconstruct=%d ",ireconstruct);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del");
   /* Read/Write axes */
   int nxv, nzv;
   float a1, a2;
   float dxv, dzv;
   float fxv, fzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("o1",a1);
   vp0.get("o2",a2);
   fxv=a2*1000.0;
   fzv=a1*1000.0;
   vp0.get("d1",a1);
   vp0.get("d2",a2);
   dzv = a1*1000.0;
   dxv = a2*1000.0;

   /* Read/Write axes from Wavefields*/
   sf_file Fx, Fz;
   sf_axis az, ax;
   int   nx, nz;
   float fx, fz;
   float dx, dz;

   Fx = sf_input("Elasticx");
   Fz = sf_input("Elasticz");

   az = sf_iaxa(Fx,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
   ax = sf_iaxa(Fx,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
   fx = sf_o(ax)*1000.0;
   fz = sf_o(az)*1000.0;

   if(nx!=nxv || nz!=nzv){
     sf_warning("Dimension not match between model and data !");
     sf_warning("nx=%d nz=%d ",nx,nz);
     sf_warning("nxv=%d nzv=%d ",nxv,nzv);
     exit(0);
   }
   if(fabs(fx-fxv)>0.1 || fabs(fz-fzv)>0.1){
     sf_warning("Coorinate original point not match between model and data !");
     sf_warning("fx=%d fz=%d ",fx,fz);
     sf_warning("fxv=%d fzv=%d ",fxv,fzv);
   }
   if(fabs(dx-dxv)>0.1 || fabs(dz-dzv)>0.1){
     sf_warning("Sampling step not match between model and data !");
     sf_warning("dx=%d dz=%d ",dx,dz);
     sf_warning("dxv=%d dzv=%d ",dxv,dzv);
   }

   sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);
   sf_warning("nx=%d nz=%d ", nx,nz);

   /* wave modeling space */
   int nxz;
   nxz=nx*nz;

   /* location of operator check */
   int ixrec=(int)((xrec1-fx)/dx);
   int izrec=(int)((zrec1-fz)/dz);
   int im1=ixrec*nz+izrec;
   sf_warning("xrec1=%f ",xrec1);
   sf_warning("zrec1=%f ",zrec1);
   sf_warning("ixrec=%d ",ixrec);
   sf_warning("izrec=%d ",izrec);
   ixrec=(int)((xrec2-fx)/dx);
   izrec=(int)((zrec2-fz)/dz);
   int im2=ixrec*nz+izrec;
   sf_warning("xrec2=%f ",xrec2);
   sf_warning("zrec2=%f ",zrec2);
   sf_warning("ixrec=%d ",ixrec);
   sf_warning("izrec=%d ",izrec);

   vp.resize(nxz);
   vs.resize(nxz);
   ep.resize(nxz);
   de.resize(nxz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;

   sf_file Fp, Fsv, Fexp1, Fexp2, Fezp1, Fezp2, Fxp1, Fxp2, Fzp1, Fzp2;

   Fp = sf_output("out");
   Fsv= sf_output("ElasticSV");

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
   float *w, *w1;
   float **ap;

   if(ireconstruct==1)
   {
   Fexp1= sf_output("Errpolxp1");
   Fexp2= sf_output("Errpolxp2");
   Fxp1= sf_output("Polxp1");
   Fxp2= sf_output("Polxp2");
   Fezp1= sf_output("Errpolzp1");
   Fezp2= sf_output("Errpolzp2");
   Fzp1= sf_output("Polzp1");
   Fzp2= sf_output("Polzp2");

   puthead2kx(Fexp1, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fexp2, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fxp1, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fxp2, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fezp1, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fezp2, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fzp1, nkz, nkx, dkz, dkx, kz0, kx0);
   puthead2kx(Fzp2, nkz, nkx, dkz, dkx, kz0, kx0);

      w = sf_floatalloc(nk);
      w1 = sf_floatalloc(nk);

      sf_warning("reconstruct x-component operators based on low-rank approx.");
      reconstruct1(w, ldataxp, fmidxp, rdataxp, nxz, nk, m2xp, n2xp, im1);

      /* error check between accurate and low-rank approx. */
      ap=sf_floatalloc2(nz, nx);

      //calculate P-wave's polarization operators directly for X-component
      polxp2dtti(ap, nx, nz, im1);

       float sum=0.0;
       int k=0;
       for(i=0;i<nx;i++)
       for(j=0;j<nz;j++)
       {
           float err=w[k]-ap[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
       }
       sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im1, sum/nxz);

       sf_floatwrite(w1, nk, Fexp1);
       sf_floatwrite(w, nk, Fxp1);

      reconstruct1(w, ldataxp, fmidxp, rdataxp, nxz, nk, m2xp, n2xp, im2);
      polxp2dtti(ap, nx, nz, im2);

       sum=0.0;
       k=0;
       for(i=0;i<nx;i++)
       for(j=0;j<nz;j++)
       {
           float err=w[k]-ap[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
       }
       sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im2, sum/nxz);

       sf_floatwrite(w1, nk, Fexp2);
       sf_floatwrite(w, nk, Fxp2);

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

   if(ireconstruct==1)
   {
      /* error check between accurate and low-rank approx. */
      reconstruct1(w, ldatazp, fmidzp, rdatazp, nxz, nk, m2zp, n2zp, im1);
      polzp2dtti(ap, nx, nz, im1);

      float sum=0.0;
      k=0;
      for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
      {
           float err=w[k]-ap[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
      }
      sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im1, sum/nxz);
      sf_floatwrite(w1, nk, Fezp1);
      sf_floatwrite(w, nk, Fzp1);

      reconstruct1(w, ldatazp, fmidzp, rdatazp, nxz, nk, m2zp, n2zp, im2);
      polzp2dtti(ap, nx, nz, im2);

      sum=0.0;
      k=0;
      for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
      {
           float err=w[k]-ap[i][j];
           sum += err*err;
           w1[k] = err;
           k++;
      }
      sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im2, sum/nxz);
      sf_floatwrite(w1, nk, Fezp2);
      sf_floatwrite(w, nk, Fzp2);

      free(*ap);
      free(w);
      free(w1);
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

   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************End of Calculating Projection Deviation Operator****************/

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

   float *px, *pz;
   float *pp, *p;

   px=sf_floatalloc(nxz);
   pz=sf_floatalloc(nxz);
   pp=sf_floatalloc(nxz);
   p=sf_floatalloc(nxz);

   /* setup I/O files */
   sf_floatread(px, nxz, Fx);
   sf_floatread(pz, nxz, Fz);

   puthead2x(Fp, nz, nx, dz/1000, dx/1000, fz/1000, fx/1000);
   puthead2x(Fsv, nz, nx, dz/1000, dx/1000, fz/1000, fx/1000);

   int iflag=0;
   sf_warning("separate qP-wave based on lowrank decomp.");
   /* separate qP wave  */
   for(k=0;k<nxz;k++) pp[k] = px[k];
   seplowrank2d(ldataxp,rdataxp,fmidxp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xp,n2xp,iflag);
   for(i=0;i<nxz;i++) p[k] = pp[k];
   for(k=0;k<nxz;k++) pp[k] = pz[k];
   seplowrank2d(ldatazp,rdatazp,fmidzp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2zp,n2zp,iflag);
   for(i=0;i<nxz;i++) p[k] += pp[k];
   sf_floatwrite(p, nxz, Fp);

   sf_warning("separate qSV-wave based on lowrank decomp.");
   /* separate qSV wave  */
   //seplowrank2d(ldataxs,rdataxs,fmidxs,px,ijkx,ijkz,nx,nz,nxz,nk,m2xs,n2xs,iflag);
   //seplowrank2d(ldatazs,rdatazs,fmidzs,pz,ijkx,ijkz,nx,nz,nxz,nk,m2zs,n2zs,iflag);
   for(i=0;i<nxz;i++) p[k] = px[k]+pz[k];
   sf_floatwrite(p, nxz, Fsv);

   t4=clock();
   timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
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

    free(px);
    free(pz);
    free(pp);
    free(p);

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

            // rotatiing according to tilted symmetry axis
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

            // rotatiing according to tilted symmetry axis
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

            // rotatiing according to tilted symmetry axis

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

            // rotatiing according to tilted symmetry axis
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

static void polxp2dtti(float **ap, int nx, int nz, int im)
{
        int    i, j, k;
        double s, c;

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
                s=sinx[k];
                c=cosx[k];
                k++;

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
        double s, c;

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
                s=sinx[k];
                c=cosx[k];
                k++;

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

