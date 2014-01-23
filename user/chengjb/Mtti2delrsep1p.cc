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
#include "ricker.h"
#include "kykxkztaper.h"
#include "fdcoef.h"
#include "eigen2x2.h"
#include "fwpttielastic.h"
#include "seplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, th;

static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

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

   float eps=0.000001;
   par.get("eps",eps,1.e-6); // tolerance
       
   int npk=20;
   par.get("npk",npk,20); // maximum rank

   int   ns=101;
   float dt=0.004, xrec=1000.0, zrec=1000.0;

   par.get("ns",ns);
   par.get("dt",dt);

   int ireconstruct=0;   // flag for reconstruct the W matrix or not
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
   sf_warning("==== Operator check point location ======");
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
    *  Calculating polarization deviation operator for wave-mode separation
    * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /********* low rank decomposition p-wave, x-component **********/
   int   m2x, n2x, m2z, n2z;
   float *ldatax, *fmidx, *rdatax;
   float *ldataz, *fmidz, *rdataz;

   iC( ddlowrank(nxz,nk,samplexp,eps,npk,lid,rid,mid) );
   m2x=mid.m();
   n2x=mid.n();
   sf_warning("m2x=%d n2x=%d",m2x, n2x);

   fmidx  = sf_floatalloc(m2x*n2x);
   ldatax = sf_floatalloc(nxz*m2x);
   rdatax = sf_floatalloc(n2x*nk);

   map2d1d(fmidx, mid, m2x, n2x);

   iC ( samplexp(md,lid,mat) );
   map2d1d(ldatax, mat, nxz, m2x);

   iC ( samplexp(rid,nd,mat) );
   map2d1d(rdatax, mat, n2x, nk);

   /********* reconsturct W-matrix for checking **********/
   int im=ixrec*nz+izrec;
   float *w;
   float **ap;
   std::valarray<float> w1(nk);
   if(ireconstruct==1)
   {
      w = sf_floatalloc(nk);

      sf_warning("reconstruct x-component operators based on low-rank approx.");
      reconstruct1(w, ldatax, fmidx, rdatax, nxz, nk, m2x, n2x, im);

      /* error check between accurate and low-rank approx. */
      ap=sf_floatalloc2(nz, nx);

      oRSF Errpolxp("Errpolxp");

      Errpolxp.put("n1",nkz);
      Errpolxp.put("n2",nkx);
      Errpolxp.put("d1",dkz);
      Errpolxp.put("d2",dkx);
      Errpolxp.put("o1",kz0);
      Errpolxp.put("o2",kx0);

      Errpolxp.put("label1","kz");
      Errpolxp.put("label2","kx");
      Errpolxp.put("unit1","2*pi/m");
      Errpolxp.put("unit2","2*pi/m");

      //calculate P-wave's polarization operators directly for X-component
      polxp2dtti(ap, nx, nz, im);

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
      sf_warning("X-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);

      Errpolxp << w1;

      oRSF Polxp("Polxp");

      Polxp.put("n1",nkz);
      Polxp.put("n2",nkx);
      Polxp.put("d1",dkz);
      Polxp.put("d2",dkx);
      Polxp.put("o1",kz0);
      Polxp.put("o2",kx0);

      Polxp.put("label1","kz");
      Polxp.put("label2","kx");
      Polxp.put("unit1","2*pi/m");
      Polxp.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++)
         w1[i]=w[i];

      Polxp << w1;

   }//ireconstruct==1

   /********* low rank decomposition p-wave, z-component **********/
   iC( ddlowrank(nxz,nk,samplezp,eps,npk,lid,rid,mid) );
   m2z=mid.m();
   n2z=mid.n();
   sf_warning("m2z=%d n2z=%d",m2z, n2z);

   fmidz  = sf_floatalloc(m2z*n2z);
   ldataz = sf_floatalloc(nxz*m2z);
   rdataz = sf_floatalloc(n2z*nk);

   map2d1d(fmidz, mid, m2z, n2z);

   iC ( samplezp(md,lid,mat) );
   map2d1d(ldataz, mat, nxz, m2z);

   iC ( samplezp(rid,nd,mat) );
   map2d1d(rdataz, mat, n2z, nk);

   /********* reconsturct W-matrix for checking **********/
   if(ireconstruct==1)
   {
      sf_warning("reconstruct z-component operators based on low-rank approx.");
      reconstruct1(w, ldataz, fmidz, rdataz, nxz, nk, m2z, n2z, im);

      /* error check between accurate and low-rank approx. */
      oRSF Errpolzp("Errpolzp");

      Errpolzp.put("n1",nkz);
      Errpolzp.put("n2",nkx);
      Errpolzp.put("d1",dkz);
      Errpolzp.put("d2",dkx);
      Errpolzp.put("o1",kz0);
      Errpolzp.put("o2",kx0);

      Errpolzp.put("label1","kz");
      Errpolzp.put("label2","kx");
      Errpolzp.put("unit1","2*pi/m");
      Errpolzp.put("unit2","2*pi/m");

      //calculate P-wave's polarization operators directly for X-component
      polzp2dtti(ap, nx, nz, im);

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
      sf_warning("Z-comp.Low-rank error: im=%d L2-norm=%20.18f",im, sum/nxz);

      Errpolzp << w1;

      oRSF Polzp("Polzp");

      Polzp.put("n1",nkz);
      Polzp.put("n2",nkx);
      Polzp.put("d1",dkz);
      Polzp.put("d2",dkx);
      Polzp.put("o1",kz0);
      Polzp.put("o2",kx0);

      Polzp.put("label1","kz");
      Polzp.put("label2","kx");
      Polzp.put("unit1","2*pi/m");
      Polzp.put("unit2","2*pi/m");

      for(i=0;i<nxz;i++)
         w1[i]=w[i];

      Polzp << w1;

      free(*ap);
      free(w);
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
   float xs, zs;
   ixs=nxv/2;
   izs=nzv/2;
   ixms=ixs+M;  /* source's x location */
   izms=izs+M;  /* source's z-location */
   xs=fx+ixs*dx;
   zs=fz+izs*dz;
   sf_warning("source location: (%f, %f)",xs, zs);

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

    float **vpp, **vss, **epp, **dee, **thee;

    vpp=sf_floatalloc2(nz,nx);
    vss=sf_floatalloc2(nz,nx);
    epp=sf_floatalloc2(nz,nx);
    dee=sf_floatalloc2(nz,nx);
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

    float *pp, *qq;
    pp=sf_floatalloc(nxz);
    qq=sf_floatalloc(nxz);

    int ii, jj, jm;

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
        fwpttielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
                      dx, dz, nx, nz, nxpad, nzpad, vpp, vss, epp, dee, thee);

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
              seplowrank2d(ldatax,rdatax,fmidx,pp,ijkx,ijkz,nx,nz,nxz,nk,m2x,n2x,iflag);
              seplowrank2d(ldataz,rdataz,fmidz,qq,ijkx,ijkz,nx,nz,nxz,nk,m2z,n2z,iflag);

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

/*********************************************************************************************************
 *            seplowrank2d(ldataxs,rdataxs,fmidxs,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xs,n2xs,iflag);
 *            seplowrank2d(ldatazs,rdatazs,fmidzs,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zs,n2zs,iflag);
 *
 *            for(i=0;i<nxz;i++)
 *              z[i]=pp[i]+qq[i];
 * *******************************************************************************************************/
              /* separate qSV wave according to that qSV's polarization vector is orthogonal to qP's  */
              sf_warning("separate qSV-wave based on lowrank decomp."); 
              seplowrank2d(ldataz,rdataz,fmidz,pp,ijkx,ijkz,nx,nz,nxz,nk,m2z,n2z,iflag);
              seplowrank2d(ldatax,rdatax,fmidx,qq,ijkx,ijkz,nx,nz,nxz,nk,m2x,n2x,iflag);

              for(i=0;i<nxz;i++)
		  z[i]=pp[i]-qq[i];   // qSV's polarization vector is orthogonal to qP's

              ElasticSepSV<<z;
          
              t5=clock();
              timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
              sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

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

    free(ldatax);
    free(ldataz);
    free(rdatax);
    free(rdataz);
    free(fmidx);
    free(fmidz);

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
    free(*thee);

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

        double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double sx = sinx[cs[b]];
            double cx = cosx[cs[b]];
            if(sx==0&&cx==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            // rotatiing according to tilted symmetry axis
            double s=sx*coss+cx*sins;
            double c=cx*coss-sx*sins;

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

        double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double sx = sinx[cs[b]];
            double cx = cosx[cs[b]];

            if(sx==0&&cx==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

            // rotatiing according to tilted symmetry axis
            double s=sx*coss+cx*sins;
            double c=cx*coss-sx*sins;

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

