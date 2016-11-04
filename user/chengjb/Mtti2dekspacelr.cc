/* 2-D two-components wavefield modeling using low-rank approximate k-space solution on the base of 
 * original elastic anisotropic displacement wave equation in TTI media.
 * But we use the original source term even if the k-space need adjustment of the
 * source.

   Copyright (C) 2014 Tongji University, Shanghai, China 
                      King Abdulah University of Science and Technology, Thuwal, Saudi Arabia
   Authors: Jiubing Cheng, Zedong Wu, and Tariq Alkhalifah
     
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
#include "eigen2x2.h"
#include "fwpvtielowrank.h"
}

static std::valarray<float> vp, vs, ep, de, th;
static double dt1, dt2, kmaxw;

static std::valarray<double>  sinx, cosx, rk2;

/* dual-domain operators based on low-rank decomp. */
int sampleux(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleuz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
int sampleuxz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

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
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");

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

   /* Fourier spectra demension */
   int nkz,nkx,nk;
   nkx=nx;
   nkz=nz;
   nk = nkx*nkz;

   /******** Stability condition given by Kosloff and Baysal (1982) for
    * Pseudo-spectral approach */
   float dttt=dx*sqrt(2.0)/(vpmax*SF_PI);
   sf_warning("Stability condition: dt < %12.8f ",dttt);

   float dkz,dkx,kz0,kx0;

   float *akx = sf_floatalloc(nk);
   float *akz = sf_floatalloc(nk);

   dkx=2*SF_PI/dx/nx;
   dkz=2*SF_PI/dz/nz;

   kx0=-SF_PI/dx;
   kz0=-SF_PI/dz;

   sf_warning("dkx=%f dkz=%f",dkx,dkz);

   sinx.resize(nk);
   cosx.resize(nk);
   rk2.resize(nk);

   double kx, kz, rk, k2;
   int    i=0, j=0, k=0, ix, iz;
   
   for(ix=0; ix < nkx; ix++)
   {
       kx = kx0+ix*dkx;

       for (iz=0; iz < nkz; iz++)
       {
            kz = kz0+iz*dkz;
            k2=kx*kx+kz*kz;
			rk=sqrt(kx*kx+kz*kz);

            rk2[i] = k2;
			akx[i] = kx;
		    akz[i] = kz;

			sinx[i] = kx/rk;
            cosx[i] = kz/rk;

            i++;
       }
   }

   /*****************************************************************************
   *  Calculating polarization deviation operator for wave-mode separation
   * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /********* low rank decomposition of operator applying to Ux **********/
   int   m2ux, n2ux, m2uz, n2uz, m2uxz, n2uxz;
   float *ldataux, *fmidux, *rdataux;
   float *ldatauz, *fmiduz, *rdatauz;
   float *ldatauxz, *fmiduxz, *rdatauxz;

   iC( ddlowrank(nxz,nk,sampleux,eps,npk,lid,rid,mid) );
   m2ux=mid.m();
   n2ux=mid.n();
   sf_warning("m2ux=%d n2ux=%d",m2ux, n2ux);

   fmidux  = sf_floatalloc(m2ux*n2ux);
   ldataux = sf_floatalloc(nxz*m2ux);
   rdataux = sf_floatalloc(n2ux*nk);

   map2d1d(fmidux, mid, m2ux, n2ux);

   iC ( sampleux(md,lid,mat) );
   map2d1d(ldataux, mat, nxz, m2ux);

   iC ( sampleux(rid,nd,mat) );
   map2d1d(rdataux, mat, n2ux, nk);
   
   /********* low rank decomposition of operator applying to Uz **********/
   iC( ddlowrank(nxz,nk,sampleuz,eps,npk,lid,rid,mid) );
   m2uz=mid.m();
   n2uz=mid.n();
   sf_warning("m2uz=%d n2uz=%d",m2uz, n2uz);

   fmiduz  = sf_floatalloc(m2uz*n2uz);
   ldatauz = sf_floatalloc(nxz*m2uz);
   rdatauz = sf_floatalloc(n2uz*nk);

   map2d1d(fmiduz, mid, m2uz, n2uz);

   iC ( sampleuz(md,lid,mat) );
   map2d1d(ldatauz, mat, nxz, m2uz);

   iC ( sampleuz(rid,nd,mat) );
   map2d1d(rdatauz, mat, n2uz, nk);

   /********* low rank decomposition of operator applying to (C13+c44)kx*kz **********/
   iC( ddlowrank(nxz,nk,sampleuxz,eps,npk,lid,rid,mid) );
   m2uxz=mid.m();
   n2uxz=mid.n();
   sf_warning("m2uxz=%d n2uxz=%d",m2uxz, n2uxz);

   fmiduxz  = sf_floatalloc(m2uxz*n2uxz);
   ldatauxz = sf_floatalloc(nxz*m2uxz);
   rdatauxz = sf_floatalloc(n2uxz*nk);

   map2d1d(fmiduxz, mid, m2uxz, n2uxz);

   iC ( sampleuxz(md,lid,mat) );
   map2d1d(ldatauxz, mat, nxz, m2uxz);

   iC ( sampleuxz(rid,nd,mat) );
   map2d1d(rdatauxz, mat, n2uxz, nk);

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
   float *p1=sf_floatalloc(nxz);
   float *p2=sf_floatalloc(nxz);
   float *p3=sf_floatalloc(nxz);
   float *pp=sf_floatalloc(nxz);

   float *q1=sf_floatalloc(nxz);
   float *q2=sf_floatalloc(nxz);
   float *q3=sf_floatalloc(nxz);
   float *qq=sf_floatalloc(nxz);

   zero1float(p1, nxz);
   zero1float(p2, nxz);
   zero1float(p3, nxz);

   zero1float(q1, nxz);
   zero1float(q2, nxz);
   zero1float(q3, nxz);

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
                  p2[(ixs+i)*nz+(izs+j)]+=i*Ricker(t, f0, t0, A);
                  q2[(ixs+i)*nz+(izs+j)]+=j*Ricker(t, f0, t0, A);
                  if(it%100==0)sf_warning("p2=%f q2=%f ",p2[(ixs+i)*nz+(izs+j)],q2[(ixs+i)*nz+(izs+j)]);
             }
        }
        /* extrapolation of Ux-componet */
   	    zero1float(pp, nxz);
        zero1float(qq, nxz);
        fwpvti2delowranksvd(ldataux,rdataux,fmidux,pp,p2,ijkx,ijkz,nx,nz,nxz,nk,m2ux,n2ux,kxm,kzm,kxzm,akx,akz);
        fwpvti2delowranksvd(ldatauxz,rdatauxz,fmiduxz,qq,q2,ijkx,ijkz,nx,nz,nxz,nk,m2uxz,n2uxz,kxm,kzm,kxzm,akx,akz);

        for(i=0;i<nxz;i++) p3[i] = pp[i] + qq[i] - p1[i];

        /******* output wavefields: components******/
        if(it==ns-1)
        {
              for(i=0;i<nxz;i++) x[i]=p3[i];
              Elasticx<<x;
        }

    	zero1float(pp, nxz);
        zero1float(qq, nxz);
        /* extrapolation of Uz-componet */
        fwpvti2delowranksvd(ldatauz,rdatauz,fmiduz,qq,q2,ijkx,ijkz,nx,nz,nxz,nk,m2uz,n2uz,kxm,kzm,kxzm,akx,akz);
        fwpvti2delowranksvd(ldatauxz,rdatauxz,fmiduxz,pp,p2,ijkx,ijkz,nx,nz,nxz,nk,m2uxz,n2uxz,kxm,kzm,kxzm,akx,akz);

        for(i=0;i<nxz;i++) q3[i] = qq[i] + pp[i] - q1[i];
        
        /******* output wavefields: components******/
        if(it==ns-1)
        {
              for(i=0;i<nxz;i++) x[i]=q3[i];
              Elasticz<<x;
        }
        /*
        if(it%100==0){
			float e;
			e = energy(p3,nxz)+energy(q3,nxz);
			sf_warning("====================Total energy: %f",e);
		}*/
        /******* update the wavefield ********/
        for(i=0;i<nxz;i++){
                p1[i]=p2[i];
                p2[i]=p3[i];

                q1[i]=q2[i];
                q2[i]=q3[i];
            }

   } //* it loop */
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for wavefield extrapolation.: %f(second)",timespent);

   free(ldataux);
   free(ldatauz);
   free(ldatauxz);
   free(rdataux);
   free(rdatauz);
   free(rdatauxz);
   free(fmidux);
   free(fmiduz);
   free(fmiduxz);

   free(p1);
   free(p2);
   free(p3);
   free(q1);
   free(q2);
   free(q3);

   free(pp);
   free(qq);

   free(ijkx);
   free(ijkz);

   free(akx);
   free(akz);
   exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* operator apply to Ux based on low-rank decomp. */
int sampleux(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44;
    double  p1, p2, s1, s2;
    double  vp2, vs2, ep2, de2, k2, kkx2, kkz2, kkxz;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;
	double  cc;
	double  a11, a22, a12;
	double  coss, sins, s, c, sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        vp2 = vp[i]*vp[i];
        vs2 = vs[i]*vs[i];
        ep2 = 1.0+2*ep[i];
		de2 = 1.0+2*de[i];
        coss=cos(th[i]);
        sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            s = sinx[cs[b]];
            c = cosx[cs[b]];
            k2=rk2[cs[b]];

            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
            }else{
              // rotatiing according to tilted symmetry axis
              sx=s*coss+c*sins;
              cx=c*coss-s*sins;
   
              kkx2 = k2*sx*sx;
              kkz2 = k2*cx*cx;

              if(sqrt(kkx2+kkz2)>=kmaxw)
                 resx(a,b) = 0.0;
              else{
                 kkxz = k2*sx*cx;

                 c33=vp2;
                 c44=vs2;
                 c11=ep2*c33;
       			 c13c44=sqrt((de2*c33-c44)*(c33-c44));

                 /* qP-wave's polarization vector */
                 a11= c11*kkx2+c44*kkz2;
                 a12= c13c44*kkxz;
                 a22= c44*kkx2+c33*kkz2;

                 aa[0][0] = a11;
                 aa[0][1] = a12;
                 aa[1][0] = a12;
                 aa[1][1] = a22;

                 dsolveSymmetric22(aa, ve, va);
 
			//sf_warning("va[0]=%f va[1]=%f ====== ok",va[0], va[1]);

                 /* qP-wave's polarization vector */
                 p1=ve[0][0];
                 p2=ve[0][1];
                 if(p1*s + p2*c <0) {
                   p1 = -p1;
                   p2 = -p2;
                 }
                 /* qS-wave's polarization vector */
                 s1=ve[1][0];
                 s2=ve[1][1];
                 if(s1*c - s2*s <0) {
                   s1 = -s1;
                   s2 = -s2;
                 }

                 lam1=sqrt(va[0])*0.5*dt1;
                 lam2=sqrt(va[1])*0.5*dt1;
                 sinclam1=sin(lam1)/lam1;
                 sinclam2=sin(lam2)/lam2;
                 sinclam1_2=sinclam1*sinclam1;
                 sinclam2_2=sinclam2*sinclam2;

				 //sf_warning("a11=%18.12f",(float)a11);
                 a11= p1*p1*va[0]*sinclam1_2 + s1*s1*va[1]*sinclam2_2;
                 //double aa11= p1*p1*va[0] + s1*s1*va[1];
				 //sf_warning("aa11=%18.12f, a11=%18.12f",(float)aa11,(float)a11);
                 cc = 2.0-dt2*a11;
            
                 resx(a,b) = cc;
             }
           }
         }// b loop
    }// a loop

    return 0;
}

/* operator apply to Uz based on low-rank decomp. */
int sampleuz(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{

    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44;
    double  p1, p2, s1, s2;
    double  vp2, vs2, ep2, de2, k2, kkx2, kkz2, kkxz;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;
	double  cc;
	double  a11, a22, a12;
	double  coss, sins, s, c, sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        vp2 = vp[i]*vp[i];
        vs2 = vs[i]*vs[i];
        ep2 = 1.0+2*ep[i];
		de2 = 1.0+2*de[i];
        coss=cos(th[i]);
        sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            s = sinx[cs[b]];
            c = cosx[cs[b]];
            k2=rk2[cs[b]];

            if(s==0&&c==0)
            {
               resx(a,b) =0.0;
            }else{
              // rotatiing according to tilted symmetry axis
              sx=s*coss+c*sins;
              cx=c*coss-s*sins;
   
              kkx2 = k2*sx*sx;
              kkz2 = k2*cx*cx;
              if(sqrt(kkx2+kkz2)>=kmaxw)
                 resx(a,b) = 0.0;
              else{
                 kkxz = k2*sx*cx;

                 c33=vp2;
                 c44=vs2;
                 c11=ep2*c33;
	     		c13c44=sqrt((de2*c33-c44)*(c33-c44));

		     	// vector decomposition operators based on polarization
                 a11= c11*kkx2+c44*kkz2;
                 a12= c13c44*kkxz;
                 a22= c44*kkx2+c33*kkz2;

                 aa[0][0] = a11;
                 aa[0][1] = a12;
                 aa[1][0] = a12;
                 aa[1][1] = a22;

                 dsolveSymmetric22(aa, ve, va);
 
                 /* qP-wave's polarization vector */
                 p1=ve[0][0];
                 p2=ve[0][1];
                 if(p1*s + p2*c <0) {
                   p1 = -p1;
                   p2 = -p2;
                 }
                 /* qS-wave's polarization vector */
                 s1=ve[1][0];
                 s2=ve[1][1];
                 if(s1*c - s2*s <0) {
                   s1 = -s1;
                   s2 = -s2;
                 }

                 lam1=sqrt(va[0])*0.5*dt1;
                 lam2=sqrt(va[1])*0.5*dt1;
                 sinclam1=sin(lam1)/lam1;
                 sinclam2=sin(lam2)/lam2;
                 sinclam1_2=sinclam1*sinclam1;
                 sinclam2_2=sinclam2*sinclam2;

				 //sf_warning("a22=%18.12f",(float)a22);
                 a22= p2*p2*va[0]*sinclam1_2 + s2*s2*va[1]*sinclam2_2;
                 //double aa22= p2*p2*va[0] + s2*s2*va[1];
				 //sf_warning("aa22=%18.12f, a22=%18.12f",(float)aa22,(float)a22);
                 cc = 2.0-dt2*a22;
            
                 resx(a,b) = cc;
               }
          }
         }// b loop
    }// a loop

    return 0;
}

int sampleuxz(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44;
    double  p1, p2, s1, s2;
    double  vp2, vs2, ep2, de2, k2, kkx2, kkz2, kkxz;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;
	double  cc;
	double  a11, a22, a12;
	double  coss, sins, s, c, sx, cx;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        vp2 = vp[i]*vp[i];
        vs2 = vs[i]*vs[i];
        ep2 = 1.0+2*ep[i];
		de2 = 1.0+2*de[i];
        coss=cos(th[i]);
        sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            s = sinx[cs[b]];
            c = cosx[cs[b]];
            k2=rk2[cs[b]];

            if(s==0&&c==0)
            {
               resx(a,b) = 0.0;
            }else{
              // rotatiing according to tilted symmetry axis
              sx=s*coss+c*sins;
              cx=c*coss-s*sins;
   
              kkx2 = k2*sx*sx;
              kkz2 = k2*cx*cx;
              if(sqrt(kkx2+kkz2)>=kmaxw)
                 resx(a,b) = 0.0;
              else{
                kkxz = k2*sx*cx;

			    // vector decomposition operators based on polarization
                c33=vp2;
                c44=vs2;
                c11=ep2*c33;
	    		c13c44=sqrt((de2*c33-c44)*(c33-c44));

                a11= c11*kkx2+c44*kkz2;
                a12= c13c44*kkxz;
                a22= c44*kkx2+c33*kkz2;

                aa[0][0] = a11;
                aa[0][1] = a12;
                aa[1][0] = a12;
                aa[1][1] = a22;

                dsolveSymmetric22(aa, ve, va);
 
                 /* qP-wave's polarization vector */
                 p1=ve[0][0];
                 p2=ve[0][1];
                 if(p1*s + p2*c <0) {
                   p1 = -p1;
                   p2 = -p2;
                 }
                 /* qS-wave's polarization vector */
                 s1=ve[1][0];
                 s2=ve[1][1];
                 if(s1*c - s2*s <0) {
                   s1 = -s1;
                   s2 = -s2;
                 }

                lam1=sqrt(va[0])*0.5*dt1;
                lam2=sqrt(va[1])*0.5*dt1;
                sinclam1=sin(lam1)/lam1;
                sinclam2=sin(lam2)/lam2;
                sinclam1_2=sinclam1*sinclam1;
                sinclam2_2=sinclam2*sinclam2;

				//sf_warning("a12=%18.12f",(float)a12);
                a12= p1*p2*va[0]*sinclam1_2 + s1*s2*va[1]*sinclam2_2;
                //double aa12= p1*p2*va[0] + s1*s2*va[1];
				//sf_warning("aa12=%18.12f, a12=%18.12f",(float)aa12,(float)a12);

                cc = -dt2*a12;
                resx(a,b) = cc;
             }
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
