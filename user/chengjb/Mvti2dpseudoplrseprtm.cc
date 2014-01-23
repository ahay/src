/* 2-D two-components wavefield modeling based on pseudo-pure mode P-wave equation
 * and P-SV separation using low-rank symbol approximation in 2D VTI media.

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
#include "puthead.h"
#include "kykxkztaper.h"
#include "fdcoef.h"
#include "eigen2x2.h"

#include "fwpvtipseudop.h"
#include "seplowrank.h"
}

/* global shared variables and arrayes */
static std::valarray<float> vp, vs, ep, de;
static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexp(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

/* theoretical deviation operator */
static void devxp2dvti(float **dev, int nx, int nz, int im);
static void devzp2dvti(float **dev, int nx, int nz, int im);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4, t5, t44;
   float   timespent, timpulse;

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
   par.get("timpulse",timpulse);

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("timpulse=%f",timpulse);
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

   sf_warning("dx=%f dz=%f fx=%f fz=%f",dx,dz,fx,fz);

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
       if(kx==0.0) kx=0.00000001*dkx;

       for (iz=0; iz < nkz; iz++)
       {
            kz = kz0+iz*dkz;
            if(kz==0.0) kz=0.00000001*dkz;
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
   int ixr, izr, ixmr, izmr;
   ixs=nxv/2;
   izs=0;
   ixms=ixs+M;  /* source's x location */
   izms=izs+M;  /* source's z-location */
   ixr=nxv/2;
   izr=0.0;
   ixmr=ixr+M;  /* receiver's x location */
   izmr=izr+M;  /* receiver's z-location */

   /* setup I/O files */
   sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6;
   Fo1 = sf_output("out");
   Fo2 = sf_output("PseudoPureSepP");
   Fo3 = sf_output("PseudoPureP1");
   Fo4 = sf_output("PseudoPureSepP1");
   Fo5 = sf_output("ImagePseudoPureP");
   Fo6 = sf_output("ImagePseudoPureSepP");

   int nss = ns/10;

   //puthead3x(Fo1, nz, nx, ns, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0f);
   //puthead3x(Fo3, nz, nx, ns, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0f);
   puthead2(Fo1, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);
   puthead2(Fo2, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);
   puthead2(Fo3, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);
   puthead2(Fo4, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);
   puthead2(Fo5, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);
   puthead2(Fo6, nz, nx, dz/1000.0, fz/1000.0, dx/1000.0, fx/1000.0);

    float *coeff_2dx=sf_floatalloc(mm);
    float *coeff_2dz=sf_floatalloc(mm);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dz,dz);

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

    float **pp1=sf_floatalloc2(nzpad, nxpad);
    float **pp2=sf_floatalloc2(nzpad, nxpad);
    float **pp3=sf_floatalloc2(nzpad, nxpad);

    float **qq1=sf_floatalloc2(nzpad, nxpad);
    float **qq2=sf_floatalloc2(nzpad, nxpad);
    float **qq3=sf_floatalloc2(nzpad, nxpad);

    zero2float(pp1, nzpad, nxpad);
    zero2float(pp2, nzpad, nxpad);
    zero2float(pp3, nzpad, nxpad);

    zero2float(qq1, nzpad, nxpad);
    zero2float(qq2, nzpad, nxpad);
    zero2float(qq3, nzpad, nxpad);

    sf_warning("==================================================");
    sf_warning("==  Porpagation Using Pseudo-Pure P-Wave Eq.    ==");
    sf_warning("==================================================");

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

    int iflag=1;

	float *image1, *image2, *z1, *z2, *zz1, *zz2;
    image1=sf_floatalloc(nxz);
    image2=sf_floatalloc(nxz);
    z1=sf_floatalloc(nxz);
    z2=sf_floatalloc(nxz);
    zz1=sf_floatalloc(nxz);
    zz2=sf_floatalloc(nxz);

    zero1float(image1, nxz);
    zero1float(image2, nxz);

    for(int it=0;it<ns;it++)
    {
	   float t=it*dt;

	   if(it%50==0)
		sf_warning("Pseudo: it= %d",it);

       /* fwpvtipseudop: forward-propagating in VTI media with pseudo-pure P-wave equation */
	   p2[ixms][izms]+=Ricker(t, f0, t0, A);
	   q2[ixms][izms]+=Ricker(t, f0, t0, A);
	   fwpvtipseudop(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz,
                      nx, nz, vpp, vss, epp, dee);

       /******* output wavefields: component and divergence *******/
       k=0;
       for(i=0;i<nx;i++)
       {
          im=i+M;
	      for(j=0;j<nz;j++)
	      {
              jm=j+M;
              pp[k] = p3[im][jm];
              qq[k] = q3[im][jm];
	          z1[k] = pp[k]+qq[k];
              k++;      
	       }
       }// i loop
	   if(it==ns-1)sf_floatwrite(z1, nxz, Fo1);

	   if(it==ns-1){
       /* correct projection error for accurate separate qP wave in spatial domain */
       sf_warning("seplowrank for x-component"); 
       seplowrank2d(ldataxp,rdataxp,fmidxp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xp,n2xp,iflag);

       sf_warning("seplowrank for z-component"); 
       seplowrank2d(ldatazp,rdatazp,fmidzp,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zp,n2zp,iflag);

       for(i=0;i<nxz;i++)
	       zz1[i]=pp[i]+qq[i];

	   sf_floatwrite(zz1, nxz, Fo2);
	   }
          
	   if(t==timpulse){
	     pp2[ixmr][0]+=1.0;
	     qq2[ixmr][0]+=1.0;
	   }
       /* bwpvtipseudop: backward-propagating in VTI media with pseudo-pure P-wave equation */
	   fwpvtipseudop(dt2, pp1, pp2, pp3, qq1, qq2, qq3, coeff_2dx, coeff_2dz,
                      nx, nz, vpp, vss, epp, dee);
	   for(ix=0;ix<nx;ix++){
	     if(t==timpulse&&ixmr==ix+M){
	        pp3[ixmr][0]=1.0;
	        qq3[ixmr][0]=1.0;
	     }else{
	        pp3[ixmr][0]=0.0;
	        qq3[ixmr][0]=0.0;
		 }
	   }

       k=0;
       for(i=0;i<nx;i++)
       {
          im=i+M;
	      for(j=0;j<nz;j++)
	      {
              jm=j+M;
              pp[k] = pp3[im][jm];
              qq[k] = qq3[im][jm];
	          z2[k] = pp[k]+qq[k];
              k++;      
	       }
       }// i loop
	   if(it==ns-1)sf_floatwrite(z2, nxz, Fo3);

       for(i=0;i<nxz;i++)
		   image1[i] += z1[i]*z2[i];

	   if(it==ns-1){
       /* correct projection error for accurate separate qP wave in spatial domain */
       sf_warning("seplowrank for x-component"); 
       seplowrank2d(ldataxp,rdataxp,fmidxp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xp,n2xp,iflag);

       sf_warning("seplowrank for z-component"); 
       seplowrank2d(ldatazp,rdatazp,fmidzp,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zp,n2zp,iflag);

       for(i=0;i<nxz;i++)
         zz2[i]=pp[i]+qq[i];

	   sf_floatwrite(zz2, nxz, Fo4);

	   }

       for(i=0;i<nxz;i++)
		   image2[i] += zz1[i]*zz2[i];

         /**************************************/
 	   for(i=0,ii=M;i<nx;i++,ii++)
	      for(j=0,jj=M;j<nz;j++,jj++)
	      {
		    p1[ii][jj]=p2[ii][jj];	
		    p2[ii][jj]=p3[ii][jj];	
		    q1[ii][jj]=q2[ii][jj];	
		    q2[ii][jj]=q3[ii][jj];	
		    pp1[ii][jj]=pp2[ii][jj];	
		    pp2[ii][jj]=pp3[ii][jj];	
		    qq1[ii][jj]=qq2[ii][jj];	
		    qq2[ii][jj]=qq3[ii][jj];	
	      }

    }/* it loop */
    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

	sf_floatwrite(image1, nxz, Fo5);
	sf_floatwrite(image2, nxz, Fo6);

    free(ldataxp);
    free(ldatazp);
    free(rdataxp);
    free(rdatazp);
    free(fmidxp);
    free(fmidzp);

    free(pp);
    free(qq);

    free(image1);
    free(image2);

    free(z1);
    free(z2);
    free(zz1);
    free(zz2);

    free(*p1);
    free(*p2);
    free(*p3);
    free(*q1);
    free(*q2);
    free(*q3);

    free(*pp1);
    free(*pp2);
    free(*pp3);
    free(*qq1);
    free(*qq2);
    free(*qq3);

    free(*vpp);
    free(*vss);
    free(*epp);
    free(*dee);

    free(coeff_2dx);
    free(coeff_2dz);

    free(ijkx);
    free(ijkz);
exit(0);
}

/* dual-domain wave-mode separation operator based on low-rank decomp. */
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
               resx(a,b) = 1.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44)); // c13+c44

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

            if(s!=0.0)
               resx(a,b) = fabs(u1/s);
            else // in fact, not happen after we move sample from zero
               resx(a,b) = 1.0;

         }// b loop
    }// a loop

    return 0;
}

/* dual-domain wave-mode separation operator based on low-rank decomp. */
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
               resz(a,b) = 1.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44)); // c13+c44

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

            if(c!=0.0)
               resz(a,b) = fabs(u2/c);
            else // in fact, not happen after we move sample from zero
               resz(a,b) = 1.0;

         }// b loop
    }// a loop

    return 0;
}
static void devxp2dvti(float **dev, int nx, int nz, int im)
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
                  dev[i][j]=1.0;

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
	       dev[i][j]=(float)(u1/sx);
          } /* j loop */
}
static void devzp2dvti(float **dev, int nx, int nz, int im)
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
                  dev[i][j]=1.0;

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
	       dev[i][j]=(float)(u2/cx);
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

