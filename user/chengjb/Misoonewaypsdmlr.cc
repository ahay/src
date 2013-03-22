/* 2-D One-way nonstationary phase-shift PSDM using low-rank approximation

   Copyright (C) 2013 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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
#include <fftw3.h>

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
/* head files aumatically produced from C programs */
extern "C"{
#include "alloc.h"
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
//#include "onewaypsdmlowrank.h"
}

static std::valarray<float> vp,wvp2;

static std::valarray<double> kx2;

static double dzz;

cpx *alloc1cpx(int n1);
cpx **alloc2cpx(int n1, int n2);
void free1cpx(cpx *p);
void free2cpx(cpx **p);

/* dual-domain phase shift one-way wave propagator based on low-rank decomp. */
int sampleup(vector<int>& rs, vector<int>& cs, CpxNumMat& resx);
int sampledo(vector<int>& rs, vector<int>& cs, CpxNumMat& resx);

static void map2d1dc(cpx *d, CpxNumMat mat, int m, int n);
static void psuplowrank2d(cpx *ldata, cpx *rdata, cpx *fmid, cpx *x, int *ijkx, int m,int n,int m2,int n2);
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
   float dt, f1, f2;

   par.get("ns",ns);
   par.get("dt",dt);
   par.get("f1",f1);
   par.get("f2",f2);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("read velocity model parameters");

   float dw=2.0*PI/(ns*dt);
   float fw=2.0*PI*f1;
   float ew=2.0*PI*f2;

   int   nw=(ew-fw)/dw+1;
   sf_warning("nw=%d dw=%f",nw,dw);
   
   /* setup I files */
   iRSF vp0;

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

   dzz=(double)dz;

   /* wave modeling space */
   int nx, nz, nxz;
   nx=nxv;
   nz=nzv;
   nxz=nx*nz;

   vp.resize(nxz);
   wvp2.resize(nx);
 
   vp0>>vp;

   /* Fourier spectra demension */
   int nk;
   nk = nx;

   float dkx,kx0;

   dkx=2*PI/dx/nx;
   kx0=-PI/dx;

   kx2.resize(nk);

   double kx, w, w2;
   int    i=0, j=0, k=0, ix, iz, ixz, iw;
   
   for(ix=0; ix < nk; ix++)
   {
       kx = kx0+ix*dkx;
       kx2[ix] = kx*kx;
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for prereparing for low-rank decomp: %f(second)",timespent);

   /* setup I/O files */
   oRSF Mig("out");
   Mig.put("n1",nz);
   Mig.put("n2",nx);
   Mig.put("d1",dz/1000);
   Mig.put("d2",dx/1000);
   Mig.put("o1",fz/1000);
   Mig.put("o2",fx/1000);

   vector<int> md(nxz), nd(nk);
   for (k=0; k < nx; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   CpxNumMat mid, mat;

   int *ijkx = sf_intalloc(nk);
   ikx(ijkx, nk);

   // set surface receiver gathers
   cpx **u, *u1;
   u = alloc2cpx(nx,nw);
   u1 = alloc1cpx(nx);

   // set image
   float **mig;
   mig = sf_floatalloc2(nz, nx);
 
   // impulse time
   float time=1.0;
   for(iw=0;iw<nw;iw++){
      w=fw+iw*dw;
      u[iw][nx/2] = cpx(0.0f, -w*time);
   }
      
   for(iz=0;iz<nz;iz++)
   {
     for(ix=0;ix<nx;ix++)
       mig[ix][iz] = 0.0;
       
     for(iw=0;iw<nw;iw++)
     {
       w = fw + iw*dw;
       w2= w*w; 
       for(ix=0;ix<nx;ix++){
         ixz=ix*nz+iz;
         //wvp2[ix] = w2/(vp[ixz]*vp[ixz]);  //prestack
         wvp2[ix] = w2/(0.25*vp[ixz]*vp[ixz]); // poststack
       }
       /********* low rank approximation of one-way phase-shift operator **********/
       int   m2up, n2up, m2do, n2do;
       cpx *ldataup, *fmidup, *rdataup;
       cpx *ldatado, *fmiddo, *rdatado;

       iC( lowrank(nx,nk,sampleup,eps,npk,lid,rid,mid) );
       m2up=mid.m();
       n2up=mid.n();
       sf_warning("m2up=%d n2up=%d",m2up, n2up);

       fmidup  = alloc1cpx(m2up*n2up);
       ldataup = alloc1cpx(nx*m2up);
       rdataup = alloc1cpx(n2up*nk);

       map2d1dc(fmidup, mid, m2up, n2up);

       iC ( sampleup(md,lid,mat) );
       map2d1dc(ldataup, mat, nxz, m2up);

       iC ( sampleup(rid,nd,mat) );
       map2d1dc(rdataup, mat, n2up, nk);

       for(ix=0;ix<nx;ix++) u1[ix]=u[iw][ix];

       psuplowrank2d(ldataup,rdataup,fmidup,u1,ijkx,nx,nk,m2up,n2up);

       // imaging condition
       for(ix=0;ix<nx;ix++){
          u[iw][ix]= u1[ix];
          mig[ix][iz] += real(u1[ix]);
       }

       free(ldataup);
       free(rdataup);
       free(fmidup);
       //free(ldatado);
       //free(rdatado);
       //free(fmiddo);
     }// iw
   }//iz
 
   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   free(ijkx);
   free2cpx(u);
   free1cpx(u1);
   free(*mig);
exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* one-way upward phase shift operator */
int sampleup(vector<int>& rs, vector<int>& cs, CpxNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,cpx(1.0f, 0.0f));

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double wv2 = wvp2[i];

        for(int b=0; b<nc; b++)
        {
            double  k2=kx2[cs[b]];

            double tmp=wv2-k2;
            if(tmp>=0.0)
               resx(a,b) = cpx(cos(-dzz*sqrt(tmp)), sin(-dzz*sqrt(tmp)));
            else
               resx(a,b) = exp(-dzz*sqrt(-tmp));
              
         }// b loop
    }// a loop

    return 0;
}

/* one-way downward phase shift operator */
int sampledo(vector<int>& rs, vector<int>& cs, CpxNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,cpx(1.0f, 0.0f));

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double wv2 = wvp2[i];

        for(int b=0; b<nc; b++)
        {
            double  k2=kx2[cs[b]];

            double tmp=wv2-k2;
            if(tmp>=0.0)
               resx(a,b) = cpx(cos(dzz*sqrt(tmp)), sin(dzz*sqrt(tmp)));
            else
               resx(a,b) = exp(-dzz*sqrt(-tmp));
              
         }// b loop
    }// a loop

    return 0;
}

static void map2d1dc(cpx *d, CpxNumMat mat, int m, int n)
{
   int i, j, k;
   k=0;
   for (i=0; i < m; i++)
   for (j=0; j < n; j++)
   {
        d[k] = mat(i,j);
        k++;
   }

}

cpx *alloc1cpx(int n1)
{
    return (cpx*)alloc1(n1,sizeof(cpx));
}

cpx **alloc2cpx(int n1, int n2)
{
        return (cpx**)alloc2(n1,n2,sizeof(cpx));
}

void free1cpx(cpx *p)
{
   free1(p);
}

void free2cpx(cpx **p)
{
   free2((void**)p);
}

/* one-way wave equation PSDM based on low-rank decomposition 
 
  Copyright (C) 2012 Tongji University (Jiubing Cheng) 
  and The University of Texas at Austin (Sergey Fomel)
 
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
/*< psuplowrank2d: one-way phase-shift migration based on low-rank decomposition >*/ 
static void psuplowrank2d(cpx *ldata, cpx *rdata, cpx *fmid, cpx *x, int *ijkx, int m,int n,int m2,int n2)
{
       int i, im, im2, jn2, ik, ikx;
       sf_complex sum1, sum2, *wp;

       wp = sf_complexalloc(m*n2);

       sf_warning("m2= %d n2=%d",m2,n2);

       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xx, *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       xx=sf_complexalloc(n);

       xp=fftwf_plan_dft_1d(m, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_1d(m,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x) to (kx) domain */
       for(i=0;i<m;i++) xin[i]=sf_cmplx(real(x[i]), imag(x[i]));

       fftwf_execute(xp);
           
       for(i=0;i<n;i++) xx[i] = xout[i];

       /* n2 IFFT from (kx) to (x) domain*/
       for(jn2=0;jn2<n2;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(ik=0;ik<m;ik++)
           {
              // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
              ikx=ijkx[ik];
              int ii=jn2n+ikx;
              xin[i]=rdata*xx;          
            }
            // (kx) to (x) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = sf_creal(xout[im])/n;
       }
       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xx);
       free(xin);
       free(xout);

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=sf_cmplx(0.0f,0.0f);
         for(im2=0;im2<m2;im2++)
         {
           sum2=sf_cmplx(0.0f,0.0f);
           for(jn2=0;jn2<n2;jn2++)
              sum2 += cmul(fmid[im2*n2+jn2],wp[jn2*m+im]);

           sum1 += cmul(ldata[im*m2+im2],sum2);
         }//im2 loop
         x[im] = sum1;
       } 

       free(wp);
}
