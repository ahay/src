/* 2-D two-components elastic wavefield extrapolation 
   using low-rank approximate PS solution on the base of 
   displacement wave equation in VTI media.

   Copyright (C) 2014 Tongji University, Shanghai, China 
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

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
#include "fwpvtielowrank.h"
#include "eigen2x2.h"
#include "decomplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, th;
static double dt1, dt2, kmaxw;

static std::valarray<double> sinx, cosx, rk2;

/* dual-domain operators based on low-rank decomp. */
int sampleax(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleaxz(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleaz(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

static void map2d1d(float *d, DblNumMat mat, int m, int n);
/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3;
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
   dt1=(double)dt;
   dt2=(double)(dt*dt);

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

   sf_warning("nx=%d nz=%d",nx,nz);
   sf_warning("dx=%f dz=%f",dx,dz);
   sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   if(nx%2==0||nz%2==0) exit(0);
   sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

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

   float kmax, vpmax, vsmax, emax;
   vpmax=0.0;
   for(int i=0;i<nxz;i++)
	   if(vp[i]>vpmax) vpmax=vp[i];

   vsmax=0.0;
   for(int i=0;i<nxz;i++)
	   if(vs[i]>vsmax) vsmax=vp[i];

   emax=0.0;
   for(int i=0;i<nxz;i++)
	   if(fabs(ep[i])>emax) emax=fabs(ep[i]);

   vpmax *= sqrt(1.0+2*emax);
   vsmax *= sqrt(1.0+2*0.3);
   sf_warning("vpmax=%f vsmax=%f emax=%f ",vpmax,vsmax,emax);

   kmax = SF_PI/dt/vsmax;
   kmaxw = (double)kmax;
   sf_warning("kmax=%f ",kmax); 

   /******** Stability condition given by Kosloff and Baysal (1982) for
	* Pseudo-spectral approach */
   float dttt=dx*sqrt(2.0)/(vpmax*SF_PI);
   sf_warning("Stability condition: dt < %12.8f ",dttt);

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

   sf_warning("dkx=%f dkz=%f",dkx,dkz);

   sinx.resize(nk);
   cosx.resize(nk);
   rk2.resize(nk);

   float *akx = sf_floatalloc(nk);
   float *akz = sf_floatalloc(nk);

   double kx, kz, rk, k2;
   int    i=0, j=0, k=0, ix, iz;
   
   for(ix=0; ix < nkx; ix++)
   {
       kx = kx0+ix*dkx;

       for (iz=0; iz < nkz; iz++)
       {
            kz = kz0+iz*dkz;

			akx[i] = kx;
            akz[i] = kz;

            k2 = kx*kx+kz*kz;
            rk = sqrt(k2);
			if(rk==0.0) rk=0.0000000001;

            sinx[i] = kx/rk;
            cosx[i] = kz/rk;
            rk2[i] = k2;
            i++;
       }
   }

   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /********* low rank decomposition of operator Ax  **********/
   int   m2ax, n2ax;
   float *ldataax, *fmidax, *rdataax;

   iC( ddlowrank(nxz,nk,sampleax,eps,npk,lid,rid,mid) );
   m2ax=mid.m();
   n2ax=mid.n();
   sf_warning("m2ax=%d n2ax=%d",m2ax, n2ax);

   fmidax  = sf_floatalloc(m2ax*n2ax);
   ldataax = sf_floatalloc(nxz*m2ax);
   rdataax = sf_floatalloc(n2ax*nk);

   map2d1d(fmidax, mid, m2ax, n2ax);

   iC ( sampleax(md,lid,mat) );
   map2d1d(ldataax, mat, nxz, m2ax);

   iC ( sampleax(rid,nd,mat) );
   map2d1d(rdataax, mat, n2ax, nk);

   /********* low rank decomposition of operator Axz  **********/
   int   m2axz, n2axz;
   float *ldataaxz, *fmidaxz, *rdataaxz;

   iC( ddlowrank(nxz,nk,sampleaxz,eps,npk,lid,rid,mid) );
   m2axz=mid.m();
   n2axz=mid.n();
   sf_warning("m2axz=%d n2axz=%d",m2axz, n2axz);

   fmidaxz  = sf_floatalloc(m2axz*n2axz);
   ldataaxz = sf_floatalloc(nxz*m2axz);
   rdataaxz = sf_floatalloc(n2axz*nk);

   map2d1d(fmidaxz, mid, m2axz, n2axz);

   iC ( sampleaxz(md,lid,mat) );
   map2d1d(ldataaxz, mat, nxz, m2axz);

   iC ( sampleaxz(rid,nd,mat) );
   map2d1d(rdataaxz, mat, n2axz, nk);

   /********* low rank decomposition of operator Az **********/
   int   m2az, n2az;
   float *ldataaz, *fmidaz, *rdataaz;

   iC( ddlowrank(nxz,nk,sampleaz,eps,npk,lid,rid,mid) );
   m2az=mid.m();
   n2az=mid.n();
   sf_warning("m2az=%d n2az=%d",m2az, n2az);

   fmidaz  = sf_floatalloc(m2az*n2az);
   ldataaz = sf_floatalloc(nxz*m2az);
   rdataaz = sf_floatalloc(n2az*nk);

   map2d1d(fmidaz, mid, m2az, n2az);

   iC ( sampleaz(md,lid,mat) );
   map2d1d(ldataaz, mat, nxz, m2az);

   iC ( sampleaz(rid,nd,mat) );
   map2d1d(rdataaz, mat, n2az, nk);

   /****************End of Calculating Projection Deviation Operator****************/
   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   /*  wavelet parameter for source definition */
   float A, f0, t0;
   f0=30.0;
   t0=0.04;
   A=1.0;

   sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);
   sf_warning("nx=%d nz=%d ", nx,nz);

   /* source definition */
   int ixs, izs;
   ixs=nxv/2;
   izs=nzv/2;

   /* setup I/O files */
   oRSF Elasticx("out");
   oRSF Elasticz("Elasticz");

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

   /********************* wavefield extrapolation *************************/
   float *ux1=sf_floatalloc(nxz);
   float *ux2=sf_floatalloc(nxz);
   float *ux3=sf_floatalloc(nxz);
   float *uz1=sf_floatalloc(nxz);
   float *uz2=sf_floatalloc(nxz);
   float *uz3=sf_floatalloc(nxz);

   float *pp=sf_floatalloc(nxz);
   float *ppp=sf_floatalloc(nxz);

   zero1float(ux1, nxz);
   zero1float(ux2, nxz);
   zero1float(ux3, nxz);
   zero1float(uz1, nxz);
   zero1float(uz2, nxz);
   zero1float(uz3, nxz);

   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikz(ijkx, ijkz, nkx, nkz);

   std::valarray<float> x(nxz);

   /* Setting Stability Conditions, by Chenlong Wang & Zedong Wu */
   float fmax = 3*f0;
   float kxm, kzm, kxzm;
   float amin, bmin;
   amin = 99999999999;
   bmin = 99999999999;
   float c11, c33, c44;
   float c1144, c3344;
    i=0;
   for (ix=0; ix<nx; ix++)
    for (iz=0; iz<nz; iz++)
    {
        c33 = vp[i] * vp[i];
        c44 = vs[i] * vs[i];
        c11 = (1+2*ep[i])*c33;
        c1144 = c11 + c44; 
        c3344 = c33 + c44;

        if (c1144<amin)
            amin = c1144;
        if (c3344<bmin)
            bmin = c3344;
        i++;
   }
   float kxmax = kx0 + nkx*dkx;
   float kzmax = kz0 + nkz*dkz;
   kxm = 2*sqrt(2)*SF_PI*fmax/sqrt(amin);
   kzm = 2*sqrt(2)*SF_PI*fmax/sqrt(bmin);
   float abmin = MIN(amin, bmin);
   kxzm = 2*sqrt(2)*SF_PI*fmax/sqrt(abmin);

   cerr<<"max kx="<<kxmax<<endl;
   cerr<<"max kz="<<kzmax<<endl;
   cerr<<"kxm="<<kxm<<endl;
   cerr<<"kzm="<<kzm<<endl;
   cerr<<"kxzm="<<kxzm<<endl;

   for(int it=0;it<ns;it++)
   {
        float t=it*dt;

        if(it%100==0)
                sf_warning("Elastic: it= %d  t=%f(s)",it,t);
 
         // 2D exploding force source
         for(i=-1;i<=1;i++)
         for(j=-1;j<=1;j++)
         {
             if(fabs(i)+fabs(j)==2)
             {
                  ux2[(ixs+i)*nz+(izs+j)]+=i*Ricker(t, f0, t0, A);
                  uz2[(ixs+i)*nz+(izs+j)]+=j*Ricker(t, f0, t0, A);
                  if(it%100==0)sf_warning("ux=%f uz=%f ",ux2[(ixs+i)*nz+(izs+j)],uz2[(ixs+i)*nz+(izs+j)]);
             }
        }
        /* extrapolation of Ux-componet */
   	    zero1float(pp, nxz);
        fwpvti2delowranksvd(ldataax,rdataax,fmidax,pp,ux2,ijkx,ijkz,nx,nz,nxz,nk,m2ax,n2ax,kxm,kzm,kxzm,akx,akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
        fwpvti2delowranksvd(ldataaxz,rdataaxz,fmidaxz,pp,uz2,ijkx,ijkz,nx,nz,nxz,nk,m2axz,n2axz,kxm,kzm,kxzm,akx,akz);
        for(i=0;i<nxz;i++)
			ux3[i] = ppp[i] - pp[i] - ux1[i];

        /* extrapolation of Uz-componet */
   	    zero1float(pp, nxz);
        fwpvti2delowranksvd(ldataaz,rdataaz,fmidaz,pp,uz2,ijkx,ijkz,nx,nz,nxz,nk,m2az,n2az,kxm,kzm,kxzm,akx,akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
        fwpvti2delowranksvd(ldataaxz,rdataaxz,fmidaxz,pp,ux2,ijkx,ijkz,nx,nz,nxz,nk,m2axz,n2axz,kxm,kzm,kxzm,akx,akz);
        for(i=0;i<nxz;i++)
			uz3[i] = ppp[i] - pp[i] - uz1[i];

        /******* update the wavefield ********/
        for(i=0;i<nxz;i++){
                ux1[i]=ux2[i];
                ux2[i]=ux3[i];
                uz1[i]=uz2[i];
                uz2[i]=uz3[i];
            }
        /******* output wavefields: components******/
        if(it==ns-1)
        {
              for(i=0;i<nxz;i++) x[i]=ux3[i];
              Elasticx<<x;
              for(i=0;i<nxz;i++) x[i]=uz3[i];
              Elasticz<<x;
        }


   } //* it loop */
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for wavefield extrapolation.: %f(second)",timespent);

   timespent=(float)(t3-t1)/(ns*CLOCKS_PER_SEC);
   sf_warning("CPU time for every time extrapolation (including low-rank decom.): %f(second)",timespent);

   free(ldataax);
   free(fmidax);
   free(rdataax);

   free(ldataaxz);
   free(fmidaxz);
   free(rdataaxz);

   free(ldataaz);
   free(fmidaz);
   free(rdataaz);

   free(ux1);
   free(ux2);
   free(ux3);
   free(uz1);
   free(uz2);
   free(uz3);

   free(pp);
   free(ppp);

   free(akx);
   free(akz);

   free(ijkx);
   free(ijkz);

   exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* operator 1 to extrapolate based on low-rank decomp. */
int sampleax(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  c44, c11, c33;
	double  sx, cx, sx2, cx2, cc;

	sf_warning("dt2=%12.9f ",dt2);

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
		double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
			// rotatiing according to tilted symmetry axis
			sx=s*coss+c*sins;
			cx=c*coss-s*sins;

            sx2=sx*sx*k2;
            cx2=cx*cx*k2;
            // wavefield extrapolator
			if(sqrt(sx2+cx2)>=kmaxw)
				resx(a,b) = 0.0;
			else{
                cc = 2.0-dt2*(c11*sx2+c44*cx2);
                resx(a,b) = cc ;
			}
              
         }// b loop
    }// a loop

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* operator 2 to extrapolate based on low-rank decomp. */
int sampleaxz(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  c44, c33, c13c44;
	double  sx, cx, kxz, cc;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double de2 = 1.0+2*de[i];
		double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
			c13c44=sqrt((de2*c33-c44)*(c33-c44));

			// rotatiing according to tilted symmetry axis
			sx=s*coss+c*sins;
			cx=c*coss-s*sins;

            // wavefield extrapolator
			if(sqrt(sx*sx*k2+cx*cx*k2)>=kmaxw)
               resx(a,b) = 0.0;
			else{
               kxz = sx*cx*k2;
               cc = c13c44*kxz*dt2;
               resx(a,b) = cc;
			}

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* operator 3 to extrapolate based on low-rank decomp. */
int sampleaz(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  c44, c33;
	double  sx, cx, sx2, cx2, cc;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
		double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;

			// rotatiing according to tilted symmetry axis
			sx=s*coss+c*sins;
			cx=c*coss-s*sins;

            sx2=sx*sx*k2;
            cx2=cx*cx*k2;

            // wavefield extrapolator
			if(sqrt(sx2+cx2)>=kmaxw)
               resx(a,b) = 0.0;
			else{
               cc = 2.0-dt2*(c44*sx2+c33*cx2);
               resx(a,b) = cc;
			}
              
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
