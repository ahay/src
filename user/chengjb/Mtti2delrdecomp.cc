/* 2-D two-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and vector decomposition based on lowrank approximation in TTI media.

   Authors: Jiubing Cheng
     
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
#include "decomplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, th;

static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation and wave vector decomposition operator based on low-rank decomp. */
static int samplexxp(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplexzp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplezzp(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplexxs(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplexzs(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplezzs(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

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

   /*****************************************************************************
   *  Calculating vector decomposition operator for P-wave
   * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   int   m2xxp, n2xxp, m2xzp, n2xzp, m2zzp, n2zzp;
   float *ldataxxp, *fmidxxp, *rdataxxp;
   float *ldataxzp, *fmidxzp, *rdataxzp;
   float *ldatazzp, *fmidzzp, *rdatazzp;

   /** Apx2 ****/
   iC( ddlowrank(nxz,nk,samplexxp,eps,npk,lid,rid,mid) );
   m2xxp=mid.m();
   n2xxp=mid.n();
   sf_warning("m2xxp=%d n2xxp=%d",m2xxp, n2xxp);

   fmidxxp  = sf_floatalloc(m2xxp*n2xxp);
   ldataxxp = sf_floatalloc(nxz*m2xxp);
   rdataxxp = sf_floatalloc(n2xxp*nk);

   map2d1d(fmidxxp, mid, m2xxp, n2xxp);

   iC ( samplexxp(md,lid,mat) );
   map2d1d(ldataxxp, mat, nxz, m2xxp);

   iC ( samplexxp(rid,nd,mat) );
   map2d1d(rdataxxp, mat, n2xxp, nk);

   /** ApxApz ****/
   iC( ddlowrank(nxz,nk,samplexzp,eps,npk,lid,rid,mid) );
   m2xzp=mid.m();
   n2xzp=mid.n();
   sf_warning("m2xzp=%d n2xzp=%d",m2xzp, n2xzp);

   fmidxzp  = sf_floatalloc(m2xzp*n2xzp);
   ldataxzp = sf_floatalloc(nxz*m2xzp);
   rdataxzp = sf_floatalloc(n2xzp*nk);

   map2d1d(fmidxzp, mid, m2xzp, n2xzp);

   iC ( samplexzp(md,lid,mat) );
   map2d1d(ldataxzp, mat, nxz, m2xzp);

   iC ( samplexzp(rid,nd,mat) );
   map2d1d(rdataxzp, mat, n2xzp, nk);

   /** Apz2 ****/
   iC( ddlowrank(nxz,nk,samplezzp,eps,npk,lid,rid,mid) );
   m2zzp=mid.m();
   n2zzp=mid.n();
   sf_warning("m2zzp=%d n2zzp=%d",m2zzp, n2zzp);

   fmidzzp  = sf_floatalloc(m2zzp*n2zzp);
   ldatazzp = sf_floatalloc(nxz*m2zzp);
   rdatazzp = sf_floatalloc(n2zzp*nk);

   map2d1d(fmidzzp, mid, m2zzp, n2zzp);

   iC ( samplezzp(md,lid,mat) );
   map2d1d(ldatazzp, mat, nxz, m2zzp);

   iC ( samplezzp(rid,nd,mat) );
   map2d1d(rdatazzp, mat, n2zzp, nk);

   /*****************************************************************************
   *  Calculating vector decomposition operator for SV-wave
   * ***************************************************************************/
   int   m2xxs, n2xxs, m2xzs, n2xzs, m2zzs, n2zzs;
   float *ldataxxs, *fmidxxs, *rdataxxs;
   float *ldataxzs, *fmidxzs, *rdataxzs;
   float *ldatazzs, *fmidzzs, *rdatazzs;

   /** Asx2 ****/
   iC( ddlowrank(nxz,nk,samplexxs,eps,npk,lid,rid,mid) );
   m2xxs=mid.m();
   n2xxs=mid.n();
   sf_warning("m2xxs=%d n2xxs=%d",m2xxs, n2xxs);

   fmidxxs  = sf_floatalloc(m2xxs*n2xxs);
   ldataxxs = sf_floatalloc(nxz*m2xxs);
   rdataxxs = sf_floatalloc(n2xxs*nk);

   map2d1d(fmidxxs, mid, m2xxs, n2xxs);

   iC ( samplexxs(md,lid,mat) );
   map2d1d(ldataxxs, mat, nxz, m2xxs);

   iC ( samplexxs(rid,nd,mat) );
   map2d1d(rdataxxs, mat, n2xxs, nk);

   /** AsxAsz ****/
   iC( ddlowrank(nxz,nk,samplexzs,eps,npk,lid,rid,mid) );
   m2xzs=mid.m();
   n2xzs=mid.n();
   sf_warning("m2xzs=%d n2xzs=%d",m2xzs, n2xzs);

   fmidxzs  = sf_floatalloc(m2xzs*n2xzs);
   ldataxzs = sf_floatalloc(nxz*m2xzs);
   rdataxzs = sf_floatalloc(n2xzs*nk);

   map2d1d(fmidxzs, mid, m2xzs, n2xzs);

   iC ( samplexzs(md,lid,mat) );
   map2d1d(ldataxzs, mat, nxz, m2xzs);

   iC ( samplexzs(rid,nd,mat) );
   map2d1d(rdataxzs, mat, n2xzs, nk);

   /** Asz2 ****/
   iC( ddlowrank(nxz,nk,samplezzs,eps,npk,lid,rid,mid) );
   m2zzs=mid.m();
   n2zzs=mid.n();
   sf_warning("m2zzs=%d n2zzs=%d",m2zzs, n2zzs);

   fmidzzs  = sf_floatalloc(m2zzs*n2zzs);
   ldatazzs = sf_floatalloc(nxz*m2zzs);
   rdatazzs = sf_floatalloc(n2zzs*nk);

   map2d1d(fmidzzs, mid, m2zzs, n2zzs);

   iC ( samplezzs(md,lid,mat) );
   map2d1d(ldatazzs, mat, nxz, m2zzs);

   iC ( samplezzs(rid,nd,mat) );
   map2d1d(rdatazzs, mat, n2zzs, nk);
   /****************End of Calculating vector Decomposition Operator****************/
   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for operator matrix low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   /*  wavelet parameter for source definition */
   float A, f0, t0;
   f0=30.0;                  
   t0=0.04;                  
   A=1.0;

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

   float xs, zs;
   xs=fx+ixs*dx;
   zs=fz+izs*dz;
   sf_warning("source location: (%f, %f)",xs, zs);

   /* setup I/O files */
   oRSF Elasticx("out"),Elasticz("ElasticZ");
   oRSF ElasticPx("ElasticPx");
   oRSF ElasticPz("ElasticPz");
   oRSF ElasticSVx("ElasticSVx");
   oRSF ElasticSVz("ElasticSVz");

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

   ElasticPx.put("n1",nkz);
   ElasticPx.put("n2",nkx);
   ElasticPx.put("d1",dz/1000);
   ElasticPx.put("d2",dx/1000);
   ElasticPx.put("o1",fz/1000);
   ElasticPx.put("o2",fx/1000);

   ElasticPz.put("n1",nkz);
   ElasticPz.put("n2",nkx);
   ElasticPz.put("d1",dz/1000);
   ElasticPz.put("d2",dx/1000);
   ElasticPz.put("o1",fz/1000);
   ElasticPz.put("o2",fx/1000);

   ElasticSVx.put("n1",nkz);
   ElasticSVx.put("n2",nkx);
   ElasticSVx.put("d1",dz/1000);
   ElasticSVx.put("d2",dx/1000);
   ElasticSVx.put("o1",fz/1000);
   ElasticSVx.put("o2",fx/1000);

   ElasticSVz.put("n1",nkz);
   ElasticSVz.put("n2",nkx);
   ElasticSVz.put("d1",dz/1000);
   ElasticSVz.put("d2",dx/1000);
   ElasticSVz.put("o1",fz/1000);
   ElasticSVz.put("o2",fx/1000);

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

    std::valarray<float> x(nxz);
    float *px, *pz;
    px=sf_floatalloc(nxz);
    pz=sf_floatalloc(nxz);

    int ii, jj;

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

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

        /* fwpttielastic: forward-propagating using original elastic equation of displacement in TTI media*/
        fwpttielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
                      dx, dz, nx, nz, nxpad, nzpad, vpp, vss, epp, dee, thee);

        /******* output wavefields: component and divergence *******/
        if(it==ns-1)
	{
              k=0;
              for(i=0;i<nx;i++){
                int im=i+M; 
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   x[k] = px[k] = p3[im][jm];
                   k++;
                }
              }
              Elasticx<<x;

              k=0;
              for(i=0;i<nx;i++){
                int im=i+M; 
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   x[k] = pz[k] = q3[im][jm];
                   k++;
                }
              }
              Elasticz<<x;

              t3=clock();
              timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
              sf_warning("CPU time for wavefield modeling: %f(second)",timespent);

              /* separate qP wave  */
              sf_warning("vector decomposition of P-wave based on lowrank decomp."); 
              decomplowrank2d(ldataxxp, rdataxxp, fmidxxp, 
                               ldataxzp, rdataxzp, fmidxzp,                              
                               ldatazzp, rdatazzp, fmidzzp,
                               px, pz, ijkx, ijkz,
                               nx, nz, nxz, nk, M, m2xxp, n2xxp, m2xzp, n2xzp, m2zzp, n2zzp);

              for(i=0;i<nxz;i++)
                 x[i] = px[i];
              ElasticPx<<x;

              for(i=0;i<nxz;i++)
                 x[i] = pz[i];
              ElasticPz<<x;
          
              k=0;
              for(i=0;i<nx;i++){
                int im=i+M; 
                for(j=0;j<nz;j++){
                   int jm=j+M;
                   px[k] = p3[im][jm];
                   pz[k] = q3[im][jm];
                   k++;
                }
              }
              /* separate qSV wave  */
              sf_warning("vector decomposition of SV-wave based on lowrank decomp."); 
              decomplowrank2d(ldataxxs, rdataxxs, fmidxxs, 
                              ldataxzs, rdataxzs, fmidxzs,                              
                              ldatazzs, rdatazzs, fmidzzs,
                              px, pz, ijkx, ijkz,
                              nx, nz, nxz, nk, M, m2xxs, n2xxs, m2xzs, n2xzs, m2zzs, n2zzs);

              for(i=0;i<nxz;i++)
                 x[i] = px[i];
              ElasticSVx<<x;

              for(i=0;i<nxz;i++)
                 x[i] = pz[i];
              ElasticSVz<<x;

              t4=clock();
              timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
              sf_warning("CPU time for vector decomposition: %f(second)",timespent);

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

    free(ldataxxp);
    free(ldatazzp);
    free(ldataxzp);

    free(rdataxxp);
    free(rdataxzp);
    free(rdatazzp);

    free(fmidxxp);
    free(fmidxzp);
    free(fmidzzp);

    free(ldataxxs);
    free(ldatazzs);
    free(ldataxzs);

    free(rdataxxs);
    free(rdataxzs);
    free(rdatazzs);

    free(fmidxxs);
    free(fmidxzs);
    free(fmidzzs);

    free(px);
    free(pz);

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

/* Apx2 */
static int samplexxp(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u1*u1;
              
         }// b loop
    }// a loop

    return 0;
}

/* ApxApz*/
static int samplexzp(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u1*u2;
              
         }// b loop
    }// a loop

    return 0;
}
/* Apz2 */ 
static int samplezzp(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            double u1=ve[0][0];
            double u2=ve[0][1];

            /* get the closest direction to k */
            if(u1*sx + u2*cx <0) {
               u1 = -u1;
               u2 = -u2;
            }

            res(a,b) = u2*u2;
              
         }// b loop
    }// a loop

    return 0;
}

/* Asx2 */
static int samplexxs(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            res(a,b) = u1*u1;
              
         }// b loop
    }// a loop

    return 0;
}

/* AsxAsz*/
static int samplexzs(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            res(a,b) = u1*u2;
              
         }// b loop
    }// a loop

    return 0;
}
/* Asz2 */ 
static int samplezzs(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

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
               res(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));

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

            res(a,b) = u2*u2;
              
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

