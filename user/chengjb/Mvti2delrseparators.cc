/* Construct low-rank approximate wave mode separators 
   for 2-D 2-C elastic anisotropic displacement wavefield
   in VTI media.

   Authors: Jiubing Cheng (Tongji University) and Sergey Fomel (The University of Texas at Austin) 
     
   Copyright (C) 2012 Tongji University(China), The University of Texas at Austin(USA) 

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

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "eigen2x2.h"
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
/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4;
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

   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);

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
   *  low-rank approximation of wave-mode separators
   * ***************************************************************************/
   /*
    *                    W = sum ( sum ( B(x,km) Amn C(xn, k) ) )
    *
    *****************************************************************************/
   /* setup I/O files */
   oRSF SepPBx("out"),   SepPAx("SepPAx"), SepPCx("SepPCx");
   oRSF SepPBz("SepPBz"),SepPAz("SepPAz"), SepPCz("SepPCz");
   oRSF SepSBx("SepSBx"),SepSAx("SepSAx"), SepSCx("SepSCx");
   oRSF SepSBz("SepSBz"),SepSAz("SepSAz"), SepSCz("SepSCz");

   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /********* low rank decomposition p-wave, x-component **********/
   int   m2, n2;
   float *lmatrix, *mmatrix, *rmatrix;
   static std::valarray<float> f1, f2, f3;

   iC( ddlowrank(nxz,nk,samplexp,eps,npk,lid,rid,mid) );
   m2=mid.m();
   n2=mid.n();
   sf_warning("X-comp of P-wave m2=%d n2=%d",m2, n2);

   mmatrix  = sf_floatalloc(m2*n2);
   lmatrix = sf_floatalloc(nxz*m2);
   rmatrix = sf_floatalloc(n2*nk);

   map2d1d(mmatrix, mid, m2, n2);

   iC ( samplexp(md,lid,mat) );
   map2d1d(lmatrix, mat, nxz, m2);

   iC ( samplexp(rid,nd,mat) );
   map2d1d(rmatrix, mat, n2, nk);
   
   f1.resize(nxz*m2);
   f2.resize(m2*n2);
   f3.resize(n2*nk);

   for(i=0;i<nxz*m2;i++)  f1[i]=lmatrix[i];
   for(i=0;i<m2*n2;i++)   f2[i]=mmatrix[i];
   for(i=0;i<n2*nk;i++)   f3[i]=rmatrix[i];

   free(lmatrix);
   free(rmatrix);
   free(mmatrix);

   SepPBx.put("n1",m2);
   SepPBx.put("n2",nxz);
   SepPBx<<f1;    // left matrix

   SepPAx.put("n1",n2);
   SepPAx.put("n2",m2);
   SepPAx<<f2;    // mid matrix

   SepPCx.put("n1",nk);
   SepPCx.put("n2",n2);
   SepPCx<<f3;    // right matrix

   /********* low rank decomposition p-wave, z-component **********/
   iC( ddlowrank(nxz,nk,samplezp,eps,npk,lid,rid,mid) );
   m2=mid.m();
   n2=mid.n();
   sf_warning("Z-comp of P-wave m2=%d n2=%d",m2, n2);

   mmatrix  = sf_floatalloc(m2*n2);
   lmatrix = sf_floatalloc(nxz*m2);
   rmatrix = sf_floatalloc(n2*nk);

   map2d1d(mmatrix, mid, m2, n2);

   iC ( samplezp(md,lid,mat) );
   map2d1d(lmatrix, mat, nxz, m2);

   iC ( samplezp(rid,nd,mat) );
   map2d1d(rmatrix, mat, n2, nk);
   
   f1.resize(nxz*m2);
   f2.resize(m2*n2);
   f3.resize(n2*nk);

   for(i=0;i<nxz*m2;i++)  f1[i]=lmatrix[i];
   for(i=0;i<m2*n2;i++)   f2[i]=mmatrix[i];
   for(i=0;i<n2*nk;i++)   f3[i]=rmatrix[i];

   free(lmatrix);
   free(rmatrix);
   free(mmatrix);

   SepPBz.put("n1",m2);
   SepPBz.put("n2",nxz);
   SepPBz<<f1;    // left matrix

   SepPAz.put("n1",n2);
   SepPAz.put("n2",m2);
   SepPAz<<f2;    // mid matrix

   SepPCz.put("n1",nk);
   SepPCz.put("n2",n2);
   SepPCz<<f3;    // right matrix

   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp for P-wave: %f(second)",timespent);

   /********* low rank decomposition SV-wave, x-component **********/
   iC( ddlowrank(nxz,nk,samplexs,eps,npk,lid,rid,mid) );
   m2=mid.m();
   n2=mid.n();
   sf_warning("X-comp of SV-wave m2=%d n2=%d",m2, n2);

   mmatrix  = sf_floatalloc(m2*n2);
   lmatrix = sf_floatalloc(nxz*m2);
   rmatrix = sf_floatalloc(n2*nk);

   map2d1d(mmatrix, mid, m2, n2);

   iC ( samplexs(md,lid,mat) );
   map2d1d(lmatrix, mat, nxz, m2);

   iC ( samplexs(rid,nd,mat) );
   map2d1d(rmatrix, mat, n2, nk);
   
   f1.resize(nxz*m2);
   f2.resize(m2*n2);
   f3.resize(n2*nk);

   for(i=0;i<nxz*m2;i++)  f1[i]=lmatrix[i];
   for(i=0;i<m2*n2;i++)   f2[i]=mmatrix[i];
   for(i=0;i<n2*nk;i++)   f3[i]=rmatrix[i];

   free(lmatrix);
   free(rmatrix);
   free(mmatrix);

   SepSBx.put("n1",m2);
   SepSBx.put("n2",nxz);
   SepSBx<<f1;    // left matrix

   SepSAx.put("n1",n2);
   SepSAx.put("n2",m2);
   SepSAx<<f2;    // mid matrix

   SepSCx.put("n1",nk);
   SepSCx.put("n2",n2);
   SepSCx<<f3;    // right matrix
   
   /********* low rank decomposition SV-wave, z-component **********/
   iC( ddlowrank(nxz,nk,samplezs,eps,npk,lid,rid,mid) );
   m2=mid.m();
   n2=mid.n();
   sf_warning("Z-comp of SV-wave m2=%d n2=%d",m2, n2);

   mmatrix  = sf_floatalloc(m2*n2);
   lmatrix = sf_floatalloc(nxz*m2);
   rmatrix = sf_floatalloc(n2*nk);

   map2d1d(mmatrix, mid, m2, n2);

   iC ( samplezs(md,lid,mat) );
   map2d1d(lmatrix, mat, nxz, m2);

   iC ( samplezs(rid,nd,mat) );
   map2d1d(rmatrix, mat, n2, nk);
   
   f1.resize(nxz*m2);
   f2.resize(m2*n2);
   f3.resize(n2*nk);

   for(i=0;i<nxz*m2;i++)  f1[i]=lmatrix[i];
   for(i=0;i<m2*n2;i++)   f2[i]=mmatrix[i];
   for(i=0;i<n2*nk;i++)   f3[i]=rmatrix[i];

   free(lmatrix);
   free(rmatrix);
   free(mmatrix);

   SepSBz.put("n1",m2);
   SepSBz.put("n2",nxz);
   SepSBz<<f1;    // left matrix

   SepSAz.put("n1",n2);
   SepSAz.put("n2",m2);
   SepSAz<<f2;    // mid matrix

   SepSCz.put("n1",nk);
   SepSCz.put("n2",n2);
   SepSCz<<f3;    // right matrix

   t4=clock();
   timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp for SV-wave: %f(second)",timespent);

   exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
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
