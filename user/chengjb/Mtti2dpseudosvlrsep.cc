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
#include "kykxkztaper.h"
#include "fdcoef.h"
#include "eigen2x2.h"

#include "fwpttipseudosv.h"
#include "seplowrank.h"
}

/* global shared variables and arrayes */
static std::valarray<float> vp, vs, ep, de, th;
static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexs(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplezs(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4, t5, t44=0;
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

   sf_warning("ns=%d dt=%f",ns,dt);
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
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

   /********* low rank decomposition p-wave, z-component **********/
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
   oRSF PseudoPureSVx("out"),PseudoPureSVz("PseudoPureSVz");
   oRSF PseudoPureSepSVx("PseudoPureSepSVx");
   oRSF PseudoPureSepSVz("PseudoPureSepSVz");
   oRSF PseudoPureSV("PseudoPureSV"),PseudoPureSepSV("PseudoPureSepSV");

   PseudoPureSVx.put("n1",nz);
   PseudoPureSVx.put("n2",nx);
   PseudoPureSVx.put("d1",dz/1000);
   PseudoPureSVx.put("d2",dx/1000);
   PseudoPureSVx.put("o1",fz/1000);
   PseudoPureSVx.put("o2",fx/1000);

   PseudoPureSVz.put("n1",nz);
   PseudoPureSVz.put("n2",nx);
   PseudoPureSVz.put("d1",dz/1000);
   PseudoPureSVz.put("d2",dx/1000);
   PseudoPureSVz.put("o1",fz/1000);
   PseudoPureSVz.put("o2",fx/1000);

   PseudoPureSV.put("n1",nz);
   PseudoPureSV.put("n2",nx);
   PseudoPureSV.put("d1",dz/1000);
   PseudoPureSV.put("d2",dx/1000);
   PseudoPureSV.put("o1",fz/1000);
   PseudoPureSV.put("o2",fx/1000);

   PseudoPureSepSV.put("n1",nz);
   PseudoPureSepSV.put("n2",nx);
   PseudoPureSepSV.put("d1",dz/1000);
   PseudoPureSepSV.put("d2",dx/1000);
   PseudoPureSepSV.put("o1",fz/1000);
   PseudoPureSepSV.put("o2",fx/1000);

   PseudoPureSepSVx.put("n1",nz);
   PseudoPureSepSVx.put("n2",nx);
   PseudoPureSepSVx.put("d1",dz/1000);
   PseudoPureSepSVx.put("d2",dx/1000);
   PseudoPureSepSVx.put("o1",fz/1000);
   PseudoPureSepSVx.put("o2",fx/1000);

   PseudoPureSepSVz.put("n1",nz);
   PseudoPureSepSVz.put("n2",nx);
   PseudoPureSepSVz.put("d1",dz/1000);
   PseudoPureSepSVz.put("d2",dx/1000);
   PseudoPureSepSVz.put("o1",fz/1000);
   PseudoPureSepSVz.put("o2",fx/1000);

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

    sf_warning("==================================================");
    sf_warning("==  Porpagation Using Pseudo-Pure P-Wave Eq.    ==");
    sf_warning("==================================================");

    std::valarray<float> x(nxz);
    std::valarray<float> y(nxz);
    std::valarray<float> z(nxz);

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

    float *pp, *qq;
    pp=sf_floatalloc(nxz);
    qq=sf_floatalloc(nxz);

    int ii, jj, im, jm;

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

    int iflag=1;

    for(int it=0;it<ns;it++)
    {
	float t=it*dt;

	if(it%50==0)
		sf_warning("Pseudo: it= %d",it);

	p2[ixms][izms]+=Ricker(t, f0, t0, A);
	q2[ixms][izms]+=Ricker(t, f0, t0, A);

        /* fwpttipseudop: forward-propagating in TTI media with pseudo-pure SV-wave equation */
	fwpttipseudosv(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz,
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

		       z[k] = x[k]+y[k];

                       k++;      
		    }
              }// i loop
              PseudoPureSVx<<x;
              PseudoPureSVz<<y;
              PseudoPureSV<<z;

              t44=clock();
              timespent=(float)(t44-t4)/CLOCKS_PER_SEC;
              sf_warning("CPU time for wavefield modeling: %f(second)",timespent);

              /* correct projection error for accurate separate qP wave in spatial domain */
              sf_warning("seplowrank for x-component"); 
              seplowrank2d(ldataxs,rdataxs,fmidxs,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xs,n2xs,iflag);

              sf_warning("seplowrank for z-component"); 
              seplowrank2d(ldatazs,rdatazs,fmidzs,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zs,n2zs,iflag);

              for(i=0;i<nxz;i++)
		  z[i]=pp[i];
              PseudoPureSepSVx<<z;
          
              for(i=0;i<nxz;i++)
		  z[i]=qq[i];
              PseudoPureSepSVz<<z;
          
              for(i=0;i<nxz;i++)
		  z[i]=pp[i]+qq[i];

              PseudoPureSepSV<<z;
          
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
    free(*thee);

    free(coeff_2dx);
    free(coeff_2dz);

    free(ijkx);
    free(ijkz);
    exit(0);
}

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexs(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

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
               resx(a,b) = 1.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44)); // c13+c44

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

            double u1=ve[1][0];
            double u2=ve[1][1];

            /* get the closest direction to k */
            if(u1*cx - u2*sx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            if(cx!=0.0)
               resx(a,b) = u1/cx;
            else // in fact, not happen after we move sample from zero
               resx(a,b) = 1.0;

         }// b loop
    }// a loop

    return 0;
}

/* dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezs(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

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
               resz(a,b) = 1.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44)); // c13+c44

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

            double u1=ve[1][0];
            double u2=ve[1][1];

            /* get the closest direction to k */
            if(u1*cx - u2*sx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            if(sx!=0.0)
               resz(a,b) = -u2/sx;
            else // in fact, not happen after we move sample from zero
               resz(a,b) = 1.0;

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

