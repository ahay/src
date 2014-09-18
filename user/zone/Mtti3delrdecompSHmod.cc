/* 3-D three-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and vector decomposition using low-rank symbol approximation in 3D TTI media.

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
#include "ricker.h"
#include "kykxkztaper.h"
#include "eigen3x3.h"
#include "puthead.h"
#include "decomplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, ga, th, ph;

static std::valarray<double> rkx, rky, rkz, rkk;

/* for low-rank decomp. */
static int samplebshxx(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebshyy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebshzz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebshxy(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebshxz(vector<int>& rs, vector<int>& cs, DblNumMat& res);
static int samplebshyz(vector<int>& rs, vector<int>& cs, DblNumMat& res);

static void map2d1d(float *d, DblNumMat mat, int m, int n);

/* definition for LAPACK SVD ROUTINEs */
char    jobz='V';  // for SVD 
char    uplo='U';  // for SVD 
int     M=3;       // for SVD 
int     LDA=M;     // for SVD 
int     LWORK=8*M; // for SVD 
int     INFO;      // for SVD 
double  Chr[9], ww[9], work[24];  // Lapack SVD array 

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

   int iflagvti;       // 1: for VTI 
   sf_warning("npk=%d ",npk);
   sf_warning("eps=%f",eps);
   sf_warning("read velocity model parameters");

   /* setup I files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), gam("gam"), the("the"), phi("phi");

   /* Read/Write axes */
   int nxv, nyv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("n3",nyv);

   float a1, a2, a3;
   vp0.get("o1",a1);
   vp0.get("o2",a2);
   vp0.get("o3",a3);

   float fxv, fyv, fzv;
   fxv=a2*1000.0;
   fyv=a3*1000.0;
   fzv=a1*1000.0;

   float dxv, dyv, dzv;
   vp0.get("d1",a1);
   vp0.get("d2",a2);
   vp0.get("d3",a3);
   dzv = a1*1000.0;
   dxv = a2*1000.0;
   dyv = a3*1000.0;

   /* Read/Write axes from Wavefields*/
   sf_file Fx, Fy, Fz;
   sf_axis az, ax, ay;
   
   Fx = sf_input("Elasticx");
   Fy = sf_input("Elasticy");
   Fz = sf_input("Elasticz");

   int   nx, ny, nz;
   float fx, fy, fz;
   float dx, dy, dz;

   az = sf_iaxa(Fx,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
   ax = sf_iaxa(Fx,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
   ay = sf_iaxa(Fx,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
   fy = sf_o(ay)*1000.0;
   fx = sf_o(ax)*1000.0;
   fz = sf_o(az)*1000.0;

   if(nx!=nxv || ny!=nyv || nz!=nzv){
     sf_warning("Dimension not match between model and data !");
     sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
     sf_warning("nxv=%d nyv=%d nzv=%d ",nxv,nyv,nzv);
     exit(0);
   }
   if(fabs(fx-fxv)>0.1 || fabs(fy-fyv)>0.1 || fabs(fz-fzv)>0.1){
     sf_warning("Coorinate original point not match between model and data !");
     sf_warning("fx=%d fy=%d fz=%d ",fx,fy,fz);
     sf_warning("fxv=%d fyv=%d fzv=%d ",fxv,fyv,fzv);
   }
   if(fabs(dx-dxv)>0.1 || fabs(dy-dyv)>0.1 || fabs(dz-dzv)>0.1){
     sf_warning("Sampling step not match between model and data !");
     sf_warning("dx=%d dy=%d dz=%d ",dx,dy,dz);
     sf_warning("dxv=%d dyv=%d dzv=%d ",dxv,dyv,dzv);
   }

   sf_warning("fx=%f fy=%f fz=%f dx=%f dy=%f dz=%f",fx,fy,fz,dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d ", nx,ny,nz);

   /* wave modeling space */
   int nxyz;
   nxyz=nx*ny*nz;

   vp.resize(nxyz);
   vs.resize(nxyz);
   ep.resize(nxyz);
   de.resize(nxyz);
   ga.resize(nxyz);
   th.resize(nxyz);
   ph.resize(nxyz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;
   gam>>ga;
   the>>th;
   phi>>ph;
  
   iflagvti=1;
   for(int i=0;i<nxyz;i++){
	  if(th[i]!=0.0||ph[i]!=0.0) iflagvti=0;
      th[i] *= SF_PI/180.0;
      ph[i] *= SF_PI/180.0;
   } 
   sf_warning("iflagvti=%d ",iflagvti);

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nx=%d ny=%d nz=%d nxyz=%d nk=%d",nx,ny,nz,nxyz,nk);

   float dkz,dkx,dky,kz0,kx0,ky0;

   dkx=2*SF_PI/dx/nx;
   dky=2*SF_PI/dy/ny;
   dkz=2*SF_PI/dz/nz;

   kx0=-SF_PI/dx;
   ky0=-SF_PI/dy;
   kz0=-SF_PI/dz;

   rkx.resize(nk);
   rky.resize(nk);
   rkz.resize(nk);
   rkk.resize(nk);

   double kx, ky, kz, k2, rk;
   int    i=0, j=0, k=0, ix, iy, iz;
   
   for(iy=0; iy < nky; iy++)
   {
     ky = ky0+iy*dky;
     //if(ky==0.0) ky=0.0000000001*dky;

     for(ix=0; ix < nkx; ix++)
     {
       kx = kx0+ix*dkx;
       //if(kx==0.0) kx=0.0000000001*dkx;

         for (iz=0; iz < nkz; iz++)
         {
            kz = kz0+iz*dkz;
            //if(kz==0.0) kz=0.0000000001*dkz;

            k2 = ky*ky+kx*kx+kz*kz;
            rk = sqrt(k2);

            rky[i] = ky/rk;
            rkx[i] = kx/rk;
            rkz[i] = kz/rk;
            rkk[i] = rk;
            i++;
         }
      }
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for prereparing for low-rank decomp: %f(second)",timespent);

   /*****************************************************************************
   *  Calculating polarization deviation operator for wave-mode separation
   * ***************************************************************************/
   vector<int> md(nxyz), nd(nk);
   for (k=0; k < nxyz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /*****************************************************************************
                    low rank approximation for vector decomposition
   * ***************************************************************************/
   int   m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, m2xy, n2xy, m2xz, n2xz, m2yz, n2yz;

   /********* qSH-bx: **********/
   iC( ddlowrank(nxyz,nk,samplebshxx,eps,npk,lid,rid,mid) );
   m2xx=mid.m();
   n2xx=mid.n();
   sf_warning("lowrank: m2xx=%d n2xx=%d",m2xx, n2xx);

   float *ldataxx, *fmidxx, *rdataxx;

   fmidxx  = sf_floatalloc(m2xx*n2xx);
   ldataxx = sf_floatalloc(nxyz*m2xx);
   rdataxx = sf_floatalloc(n2xx*nk);

   map2d1d(fmidxx, mid, m2xx, n2xx);

   iC ( samplebshxx(md,lid,mat) );
   map2d1d(ldataxx, mat, nxyz, m2xx);

   iC ( samplebshxx(rid,nd,mat) );
   map2d1d(rdataxx, mat, n2xx, nk);

   /********* qSH-by: **********/
   iC( ddlowrank(nxyz,nk,samplebshyy,eps,npk,lid,rid,mid) );
   m2yy=mid.m();
   n2yy=mid.n();
   sf_warning("lowrank: m2yy=%d n2yy=%d",m2yy, n2yy);

   float *ldatayy, *fmidyy, *rdatayy;

   fmidyy  = sf_floatalloc(m2yy*n2yy);
   ldatayy = sf_floatalloc(nxyz*m2yy);
   rdatayy = sf_floatalloc(n2yy*nk);

   map2d1d(fmidyy, mid, m2yy, n2yy);

   iC ( samplebshyy(md,lid,mat) );
   map2d1d(ldatayy, mat, nxyz, m2yy);

   iC ( samplebshyy(rid,nd,mat) );
   map2d1d(rdatayy, mat, n2yy, nk);

   float *ldatazz, *fmidzz, *rdatazz;
   /********* qSH-bz: **********/
   if(iflagvti!=1){
   iC( ddlowrank(nxyz,nk,samplebshzz,eps,npk,lid,rid,mid) );
   m2zz=mid.m();
   n2zz=mid.n();
   sf_warning("lowrank: m2zz=%d n2zz=%d",m2zz, n2zz);

   fmidzz  = sf_floatalloc(m2zz*n2zz);
   ldatazz = sf_floatalloc(nxyz*m2zz);
   rdatazz = sf_floatalloc(n2zz*nk);

   map2d1d(fmidzz, mid, m2zz, n2zz);

   iC ( samplebshzz(md,lid,mat) );
   map2d1d(ldatazz, mat, nxyz, m2zz);

   iC ( samplebshzz(rid,nd,mat) );
   map2d1d(rdatazz, mat, n2zz, nk);
   }

   /********* AxAy-operator **********/
   iC( ddlowrank(nxyz,nk,samplebshxy,eps,npk,lid,rid,mid) );
   m2xy=mid.m();
   n2xy=mid.n();
   sf_warning("lowrank: m2xy=%d n2xy=%d",m2xy, n2xy);

   float *ldataxy, *fmidxy, *rdataxy;

   fmidxy  = sf_floatalloc(m2xy*n2xy);
   ldataxy = sf_floatalloc(nxyz*m2xy);
   rdataxy = sf_floatalloc(n2xy*nk);

   map2d1d(fmidxy, mid, m2xy, n2xy);

   iC ( samplebshxy(md,lid,mat) );
   map2d1d(ldataxy, mat, nxyz, m2xy);

   iC ( samplebshxy(rid,nd,mat) );
   map2d1d(rdataxy, mat, n2xy, nk);

   float *ldataxz, *fmidxz, *rdataxz;
   float *ldatayz, *fmidyz, *rdatayz;

   if(iflagvti!=1){
   /********* AxAz-operator **********/
   iC( ddlowrank(nxyz,nk,samplebshxz,eps,npk,lid,rid,mid) );
   m2xz=mid.m();
   n2xz=mid.n();
   sf_warning("lowrank: m2xz=%d n2xz=%d",m2xz, n2xz);

   fmidxz  = sf_floatalloc(m2xz*n2xz);
   ldataxz = sf_floatalloc(nxyz*m2xz);
   rdataxz = sf_floatalloc(n2xz*nk);

   map2d1d(fmidxz, mid, m2xz, n2xz);

   iC ( samplebshxz(md,lid,mat) );
   map2d1d(ldataxz, mat, nxyz, m2xz);

   iC ( samplebshxz(rid,nd,mat) );
   map2d1d(rdataxz, mat, n2xz, nk);

   /********* AyAz-operator **********/
   iC( ddlowrank(nxyz,nk,samplebshyz,eps,npk,lid,rid,mid) );
   m2yz=mid.m();
   n2yz=mid.n();
   sf_warning("lowrank: m2yz=%d n2yz=%d",m2yz, n2yz);

   fmidyz  = sf_floatalloc(m2yz*n2yz);
   ldatayz = sf_floatalloc(nxyz*m2yz);
   rdatayz = sf_floatalloc(n2yz*nk);

   map2d1d(fmidyz, mid, m2yz, n2yz);

   iC ( samplebshyz(md,lid,mat) );
   map2d1d(ldatayz, mat, nxyz, m2yz);

   iC ( samplebshyz(rid,nd,mat) );
   map2d1d(rdatayz, mat, n2yz, nk);

   }

   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to calculate wavefield****************/
   /****************begin to calculate wavefield****************/
   //  wavelet parameter for source definition */

   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

    float *px, *py, *pz, *x;
    px=sf_floatalloc(nxyz);
    py=sf_floatalloc(nxyz);
    pz=sf_floatalloc(nxyz);
    x=sf_floatalloc(nxyz);

    int iflag;

    sf_floatread(px, nxyz, Fx);
    sf_floatread(py, nxyz, Fy);
    sf_floatread(pz, nxyz, Fz);

    sf_file Fsx, Fsy, Fsz;
    Fsx= sf_output("out");
    Fsy= sf_output("ElasticSHy");
    Fsz= sf_output("ElasticSHz");

    puthead3x(Fsx, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
    puthead3x(Fsy, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
    puthead3x(Fsz, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);

    // separate and decompose qSH wave 
    sf_warning("vector decomposition of SH-wave based on lowrank decomp.");
	if(iflagvti!=1)
       decomplowrank3dp(ldataxx, rdataxx, fmidxx,
                    ldatayy, rdatayy, fmidyy,
                    ldatazz, rdatazz, fmidzz,
                    ldataxy, rdataxy, fmidxy,
                    ldataxz, rdataxz, fmidxz,
                    ldatayz, rdatayz, fmidyz,
                    px, py, pz, ijkx, ijky, ijkz,
                    nx, ny, nz, nxyz, nk,
                    m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, 
                    m2xz, n2xz, m2xy, n2xy, m2yz, n2yz);
	else{
       decomplowrank3dshvti(ldataxx, rdataxx, fmidxx,
                    ldatayy, rdatayy, fmidyy,
                    ldataxy, rdataxy, fmidxy,
                    px, py, ijkx, ijky, ijkz,
                    nx, ny, nz, nxyz, nk,
                    m2xx, n2xx, m2yy, n2yy, m2zz, n2zz, 
                    m2xy, n2xy);
	   zero1float(pz, nxyz);
	}


	sf_floatwrite(px, nxyz, Fsx);
	sf_floatwrite(py, nxyz, Fsy);
	sf_floatwrite(pz, nxyz, Fsz);

    free(px);
    free(py);
    free(pz);

    t5=clock();
    timespent=(float)(t5-t44)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes decomposition.: %f(second)",timespent);

    free(ldataxx);
    free(ldatayy);
    free(rdataxx);
    free(rdatayy);
    free(fmidxx);
    free(fmidyy);

	if(iflagvti!=1)
	{
      free(ldatazz);
      free(rdatazz);
      free(fmidzz);
      free(ldataxy);
      free(ldataxz);
      free(ldatayz);
      free(rdataxy);
      free(rdataxz);
      free(rdatayz);
      free(fmidxy);
      free(fmidxz);
      free(fmidyz);
    }

    sf_warning("-------sucessful ending --------");
    exit(0);
}

/* Asvy2  */
static int samplebshyy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss;
            cc = kx*xn+ky*yn+kz*zn;
            ss = 1-cc*cc;

            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usy /= rk;
               res(a,b) = usy*usy*ss;
            }

         }// b loop
    }// a loop

    return 0;
}

/* Asvx2 */
static int samplebshxx(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss;
            cc = kx*xn+ky*yn+kz*zn;
            ss = 1-cc*cc;

            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usx /= rk;
               res(a,b) = usx*usx*ss;
            }

         }// b loop
    }// a loop

    return 0;
}

/* Asvz2 */
static int samplebshzz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss;
            cc = kx*xn+ky*yn+kz*zn;
            ss = 1-cc*cc;

            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usz /= rk;
               res(a,b) = usz*usz*ss;
            }
              
         }// b loop
    }// a loop

    return 0;
}

/* AsvyAsvx  */
static int samplebshxy(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);
 
    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss;
            cc = kx*xn+ky*yn+kz*zn;
            ss = sqrt(1-cc*cc);

            double rk=sqrt(usx*usx+usy*usy+usz*usz);

            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usx /= rk;
               usy /= rk;
               res(a,b) = usx*usy*ss;
            }

         }// b loop
    }// a loop

    return 0;
}

/* AsvxAsvz  */
static int samplebshxz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss; 
            cc = kx*xn+ky*yn+kz*zn;
            ss = 1-cc*cc;

            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usx /= rk;
               usz /= rk;
               res(a,b) = usx*usz*ss;
            }

         }// b loop
    }// a loop

    return 0;
}

/* AsvyAsvz  */
static int samplebshyz(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();

    res.resize(nr,nc);

    setvalue(res,0.0);

    double usx, usy, usz;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        double cos_ph=cos(ph[i]);
        double sin_ph=sin(ph[i]);
        double cos_th=cos(th[i]);
        double sin_th=sin(th[i]);

        // define axis vector
        double xn= cos_ph*sin_th;
        double yn= sin_ph*sin_th;
        double zn= cos_th;

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               res(a,b) = 0.0;
               continue;
            }

            /* define SH's polarization in TTI medium */
            usx= kz*yn - ky*zn;
            usy= kx*zn - kz*xn;
            usz= ky*xn - kx*yn;

            /* define the polar angle*/
            double cc, ss;
            cc = kx*xn+ky*yn+kz*zn;
            ss = sqrt(1-cc*cc);

            double rk=sqrt(usx*usx+usy*usy+usz*usz);
            if(rk==0.0)
               res(a,b) = 0.0;
            else{
               usy /= rk;
               usz /= rk;
               res(a,b) = usy*usz*ss;
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
