/* vector decomposition of the inputed elastic wavefields based on lowrank approximation in TTI media 

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
#include <rsf.h>
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
#include "kykxkztaper.h"
#include "eigen2x2.h"
#include "decomplowrank.h"
#include "seplowrank.h"  // for sub-program: reconstruct(...)
}

static std::valarray<float> vp, vs, ep, de, th;

static std::valarray<double> sinx, cosx, rkk;

/* dual-domain wave-mode separation and wave vector decomposition operator based on low-rank decomp. */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

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
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), the("the");

   /* Read/Write axes */
   int nxv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);

   float a1, a2;
   vp0.get("o1",a1);
   vp0.get("o2",a2);

   float fxv, fzv;
   fxv=a2*1000.0;
   fzv=a1*1000.0;

   float dxv, dzv;
   vp0.get("d1",a1);
   vp0.get("d2",a2);
   dxv = a1*1000.0;
   dzv = a2*1000.0;

   /* Read/Write axes from Wavefields*/
   sf_file Fx, Fz;
   sf_axis az, ax;
   
   Fx = sf_input("Elasticx");
   Fz = sf_input("Elasticz");

   int   nx, nz;
   float fx, fz;
   float dx, dz;
   
   az = sf_iaxa(Fx,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
   ax = sf_iaxa(Fx,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
   fx = sf_o(ax)*1000.0;
   fz = sf_o(az)*1000.0;

   sf_warning("nx=%d nz=%d ",nx,nz);
   sf_warning("nxv=%d nzv=%d ",nxv,nzv);
   sf_warning("fx=%f fz=%f ",fx,fz);
   sf_warning("fxv=%f fzv=%f ",fxv,fzv);
   sf_warning("dx=%f dz=%f ",dx,dz);
   sf_warning("dxv=%f dzv=%f ",dxv,dzv);

   if(nx!=nxv || nz!=nzv){
     sf_warning("Dimension not match between model and data !");
     exit(0);
   }
   if(fabs(fx-fxv)>0.1 || fabs(fz-fzv)>0.1){
     sf_warning("Coorinate original point not match between model and data !");
   }
   if(fabs(dx-dxv)>0.1 || fabs(dz-dzv)>0.1){
     sf_warning("Sampling step not match between model and data !");
   }
   
   /* wave modeling space */
   int nxz;
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
   *  Calculating vector decomposition operator for P-wave
   * ***************************************************************************/
   vector<int> md(nxz), nd(nk);
   for (k=0; k < nxz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   int   m2xx, n2xx, m2xz, n2xz, m2zz, n2zz;
   float *ldataxx, *fmidxx, *rdataxx;
   float *ldataxz, *fmidxz, *rdataxz;
   float *ldatazz, *fmidzz, *rdatazz;

   /** Apx2 ****/
   iC( ddlowrank(nxz,nk,samplexx,eps,npk,lid,rid,mid) );
   m2xx=mid.m();
   n2xx=mid.n();
   sf_warning("m2xx=%d n2xx=%d",m2xx, n2xx);

   fmidxx  = sf_floatalloc(m2xx*n2xx);
   ldataxx = sf_floatalloc(nxz*m2xx);
   rdataxx = sf_floatalloc(n2xx*nk);

   map2d1d(fmidxx, mid, m2xx, n2xx);

   iC ( samplexx(md,lid,mat) );
   map2d1d(ldataxx, mat, nxz, m2xx);

   iC ( samplexx(rid,nd,mat) );
   map2d1d(rdataxx, mat, n2xx, nk);

  /** ApxApz ****/
   iC( ddlowrank(nxz,nk,samplexz,eps,npk,lid,rid,mid) );
   m2xz=mid.m();
   n2xz=mid.n();
   sf_warning("m2xz=%d n2xz=%d",m2xz, n2xz);

   fmidxz  = sf_floatalloc(m2xz*n2xz);
   ldataxz = sf_floatalloc(nxz*m2xz);
   rdataxz = sf_floatalloc(n2xz*nk);

   map2d1d(fmidxz, mid, m2xz, n2xz);

   iC ( samplexz(md,lid,mat) );
   map2d1d(ldataxz, mat, nxz, m2xz);

   iC ( samplexz(rid,nd,mat) );
   map2d1d(rdataxz, mat, n2xz, nk);

   /** Apz2 ****/
   iC( ddlowrank(nxz,nk,samplezz,eps,npk,lid,rid,mid) );
   m2zz=mid.m();
   n2zz=mid.n();
   sf_warning("m2zz=%d n2zz=%d",m2zz, n2zz);

   fmidzz  = sf_floatalloc(m2zz*n2zz);
   ldatazz = sf_floatalloc(nxz*m2zz);
   rdatazz = sf_floatalloc(n2zz*nk);

   map2d1d(fmidzz, mid, m2zz, n2zz);

   iC ( samplezz(md,lid,mat) );
   map2d1d(ldatazz, mat, nxz, m2zz);

   iC ( samplezz(rid,nd,mat) );
   map2d1d(rdatazz, mat, n2zz, nk);

   /****************End of Calculating vector Decomposition Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   oRSF ElasticPx("ElasticPx");
   oRSF ElasticPz("ElasticPz");
   oRSF ElasticSx("ElasticSx");
   oRSF ElasticSz("ElasticSz");
   oRSF Orthog("Orthog");

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

   ElasticSx.put("n1",nkz);
   ElasticSx.put("n2",nkx);
   ElasticSx.put("d1",dz/1000);
   ElasticSx.put("d2",dx/1000);
   ElasticSx.put("o1",fz/1000);
   ElasticSx.put("o2",fx/1000);

   ElasticSz.put("n1",nkz);
   ElasticSz.put("n2",nkx);
   ElasticSz.put("d1",dz/1000);
   ElasticSz.put("d2",dx/1000);
   ElasticSz.put("o1",fz/1000);
   ElasticSz.put("o2",fx/1000);

   Orthog.put("n1",nkz);
   Orthog.put("n2",nkx);
   Orthog.put("d1",dz/1000);
   Orthog.put("d2",dx/1000);
   Orthog.put("o1",fz/1000);
   Orthog.put("o2",fx/1000);

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
    float *px, *pz, *px1, *pz1;
    px=sf_floatalloc(nxz);
    pz=sf_floatalloc(nxz);
    px1=sf_floatalloc(nxz);
    pz1=sf_floatalloc(nxz);

    sf_floatread(px, nxz, Fx);
    sf_floatread(pz, nxz, Fz);

    for(i=0;i<nxz;i++) {
		px1[i] = px[i]; 
		pz1[i] = pz[i]; 
	}

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

	int M=0;

    /* separate qP wave  */
    sf_warning("vector decomposition of P-wave based on lowrank decomp."); 
    decomplowrank2d(ldataxx, rdataxx, fmidxx, 
                             ldataxz, rdataxz, fmidxz,                              
                             ldatazz, rdatazz, fmidzz,
                             px, pz, ijkx, ijkz,
                             nx, nz, nxz, nk, M, m2xx, n2xx, m2xz, n2xz, m2zz, n2zz);

    sf_warning("vector decomposition of SV-wave based on lowrank decomp."); 
    decomplowrank2ds(ldataxx, rdataxx, fmidxx, 
                              ldataxz, rdataxz, fmidxz,                              
                              ldatazz, rdatazz, fmidzz,
                              px1, pz1, ijkx, ijkz,
                              nx, nz, nxz, nk, M, m2xx, n2xx, m2xz, n2xz, m2zz, n2zz);

    t4=clock();
    timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

    for(i=0;i<nxz;i++) x[i] = px[i]; 
    ElasticPx<<x;

    for(i=0;i<nxz;i++) x[i] = pz[i];
    ElasticPz<<x;

    for(i=0;i<nxz;i++) x[i] = px1[i];
    ElasticSx<<x;

    for(i=0;i<nxz;i++) x[i] = pz1[i];
    ElasticSz<<x;

    for(i=0;i<nxz;i++) x[i] = px[i]*px1[i]+pz[i]*pz1[i];
    Orthog<<x;

    free(ldataxx);
    free(ldatazz);
    free(ldataxz);

    free(rdataxx);
    free(rdataxz);
    free(rdatazz);

    free(fmidxx);
    free(fmidxz);
    free(fmidzz);

    free(px);
    free(pz);
    free(px1);
    free(pz1);

    free(*vpp);
    free(*vss);
    free(*epp);
    free(*dee);
    free(*thee);

    free(ijkx);
    free(ijkz);

    exit(0);
}

/* Apx2 */
static int samplexx(vector<int>& rs, vector<int>& cs, DblNumMat& res)
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
static int samplexz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
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
static int samplezz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
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
