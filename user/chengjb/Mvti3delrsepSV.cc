/* 3-D three-components wavefield modeling based on original elastic anisotropic displacement 
  wave equation and P-S separation using low-rank symbol approximation in 3D VTI media.

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
#include "seplowrank.h"
#include "puthead.h"
}

static std::valarray<float> vp, vs, ep, de, ga;

static std::valarray<double> rkx, rky, rkz;

/* for low-rank decomp. */

static int samplexsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
static int sampleysv3(vector<int>& rs, vector<int>& cs, DblNumMat& resy);
static int samplezsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resz);

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

   /* setup I/O files */
   iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), gam("gam");

   float a1, a2, a3;

   /* Read/Write axes from Model*/
   int nxv, nyv, nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("n3",nyv);
   vp0.get("o1",a1);
   vp0.get("o2",a2);
   vp0.get("o3",a3);
   float fxv, fyv, fzv;
   fxv=a1*1000.0;
   fyv=a2*1000.0;
   fzv=a3*1000.0;
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
   int nxyz=nx*ny*nz;

   vp.resize(nxyz);
   vs.resize(nxyz);
   ep.resize(nxyz);
   de.resize(nxyz);
   ga.resize(nxyz);
 
   vp0>>vp;
   vs0>>vs;
   epsi>>ep;
   del>>de;
   gam>>ga;

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nxyz=%d nk=%d",nxyz,nk);

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

   double kx, ky, kz, k2, rk;
   int    i=0, k=0, ix, iy, iz;
   
   for(iy=0; iy < nky; iy++)
   {
     ky = ky0+iy*dky;
     if(ky==0.0) ky=0.0000000001*dky;

     for(ix=0; ix < nkx; ix++)
     {
       kx = kx0+ix*dkx;
       if(kx==0.0) kx=0.0000000001*dkx;

         for (iz=0; iz < nkz; iz++)
         {
            kz = kz0+iz*dkz;
            if(kz==0.0) kz=0.0000000001*dkz;

            k2 = ky*ky+kx*kx+kz*kz;
            rk = sqrt(k2);

            rky[i] = ky/rk;
            rkx[i] = kx/rk;
            rkz[i] = kz/rk;
            i++;
         }
      }
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;

   /*****************************************************************************
   *  Calculating polarization deviation operator for wave-mode separation
   * ***************************************************************************/
   vector<int> md(nxyz), nd(nk);
   for (k=0; k < nxyz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   /*****************************************************************************
                    low rank decomposition for SV-wave's operators
   * ***************************************************************************/
   int   m2ysv, n2ysv, m2xsv, n2xsv, m2zsv, n2zsv;
   float *ldataxsv, *fmidxsv, *rdataxsv;
   float *ldataysv, *fmidysv, *rdataysv;
   float *ldatazsv, *fmidzsv, *rdatazsv;

   /********* low rank decomposition SV-wave, y-component **********/
   iC( ddlowrank(nxyz,nk,sampleysv3,eps,npk,lid,rid,mid) );
   m2ysv=mid.m();
   n2ysv=mid.n();
   sf_warning("lowrank-sv-y:m2ysv=%d n2ysv=%d",m2ysv, n2ysv);

   fmidysv  = sf_floatalloc(m2ysv*n2ysv);
   ldataysv = sf_floatalloc(nxyz*m2ysv);
   rdataysv = sf_floatalloc(n2ysv*nk);

   map2d1d(fmidysv, mid, m2ysv, n2ysv);

   iC ( sampleysv3(md,lid,mat) );
   map2d1d(ldataysv, mat, nxyz, m2ysv);

   iC ( sampleysv3(rid,nd,mat) );
   map2d1d(rdataysv, mat, n2ysv, nk);

   /********* low rank decomposition SV-wave, x-component **********/

   iC( ddlowrank(nxyz,nk,samplexsv3,eps,npk,lid,rid,mid) );
   m2xsv=mid.m();
   n2xsv=mid.n();
   sf_warning("lowrank-sv-x:m2xsv=%d n2xsv=%d",m2xsv, n2xsv);

   fmidxsv  = sf_floatalloc(m2xsv*n2xsv);
   ldataxsv = sf_floatalloc(nxyz*m2xsv);
   rdataxsv = sf_floatalloc(n2xsv*nk);

   map2d1d(fmidxsv, mid, m2xsv, n2xsv);

   iC ( samplexsv3(md,lid,mat) );
   map2d1d(ldataxsv, mat, nxyz, m2xsv);

   iC ( samplexsv3(rid,nd,mat) );
   map2d1d(rdataxsv, mat, n2xsv, nk);

   /********* low rank decomposition SV-wave, z-component **********/

   iC( ddlowrank(nxyz,nk,samplezsv3,eps,npk,lid,rid,mid) );
   m2zsv=mid.m();
   n2zsv=mid.n();
   sf_warning("lowrank-sv-z:m2zsv=%d n2zsv=%d",m2zsv, n2zsv);

   fmidzsv  = sf_floatalloc(m2zsv*n2zsv);
   ldatazsv = sf_floatalloc(nxyz*m2zsv);
   rdatazsv = sf_floatalloc(n2zsv*nk);

   map2d1d(fmidzsv, mid, m2zsv, n2zsv);

   iC ( samplezsv3(md,lid,mat) );
   map2d1d(ldatazsv, mat, nxyz, m2zsv);

   iC ( samplezsv3(rid,nd,mat) );
   map2d1d(rdatazsv, mat, n2zsv, nk);
   
   /****************End of Calculating Projection Deviation Operator****************/
   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

   /****************begin to input wavefield (snapshot) and separation ****************/
   /****************begin to input wavefield (snapshot) and separation ****************/
   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

    float *pp, *p;

    pp=sf_floatalloc(nxyz);
    p=sf_floatalloc(nxyz);

    int iflag;

	sf_file F;
    F = sf_output("out");
	puthead3x(F, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);

    sf_warning("separate qSV-wave based on lowrank decomp."); 
    // separate qSV wave  
    iflag=1;   /*because Asv is a cross product of Ap and Ash, there has a i */
    for(k=0;k<nxyz;k++) p[k] = 0.0;

	sf_floatread(pp, nxyz, Fx);
    seplowrank3d(ldataxsv,rdataxsv,fmidxsv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2xsv,n2xsv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

	sf_floatread(pp, nxyz, Fy);
    seplowrank3d(ldataysv,rdataysv,fmidysv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2ysv,n2ysv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

	sf_floatread(pp, nxyz, Fz);
    seplowrank3d(ldatazsv,rdatazsv,fmidzsv,pp,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2zsv,n2zsv,iflag);
    for(k=0;k<nxyz;k++) p[k] += pp[k];

	sf_floatwrite(p, nxyz, F);

    free(pp);
    free(p);

    t4=clock();
    timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

    free(ldataxsv);
    free(ldataysv);
    free(ldatazsv);
    free(rdataxsv);
    free(rdataysv);
    free(rdatazsv);
    free(fmidxsv);
    free(fmidysv);
    free(fmidzsv);
    sf_warning("-------sucessful ending --------");
    exit(0);
}

/* SV-wave y-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int sampleysv3(vector<int>& rs, vector<int>& cs, DblNumMat& resy)
{
    int nr = rs.size();
    int nc = cs.size();

    resy.resize(nr,nc);

    setvalue(resy,0.0);

    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               resy(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
                upy=-Chr[7];
                upz=-Chr[8];
            }

            double usx= -kx*upz;
            double usy= -ky*upz;
            double usz= kx*upx + ky*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
			if(rk==0.0)
               resy(a,b) = 0.0;
			else{
               usy /= rk;
			   double ss=sqrt(1.0-kz*kz);
               resy(a,b) = usy*ss;   /* scaling operator to tackle kiss singularity */
			}
              
         }// b loop
    }// a loop

    return 0;
}

/* SV-wave x-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplexsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               resx(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
                upy=-Chr[7];
                upz=-Chr[8];
            }

            double usx= -kx*upz;
            double usy= -ky*upz;
            double usz= kx*upx + ky*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
			if(rk==0.0)
                resx(a,b) = 0.0;
			else{
                usx /= rk;
			    double ss=sqrt(1.0-kz*kz);
                resx(a,b) = usx*ss;   /* scaling operator to tackle kiss singularity */
			}
              
         }// b loop
    }// a loop

    return 0;
}
/* SV-wave z-component dual-domain wave-mode separation operator based on low-rank decomp. */
static int samplezsv3(vector<int>& rs, vector<int>& cs, DblNumMat& resz)
{
    int nr = rs.size();
    int nc = cs.size();

    resz.resize(nr,nc);

    setvalue(resz,0.0);

    double upx, upy, upz;

    double c33, c44, c66, c11, c13c44, c11c66, a11, a12, a22, a33, a13, a23;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
        double ga2 = 1.0+2*ga[i];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            if(kx==0&&ky==0&&kz==0)
            {
               resz(a,b) = 0.0;
               continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
            c66=ga2*c44;
            c13c44=sqrt((de2*c33-c44)*(c33-c44));
            c11c66=c11-c66;

            a11= c11*kx*kx + c66*ky*ky + c44*kz*kz;
            a22= c66*kx*kx + c11*ky*ky + c44*kz*kz;
            a33= c44*kx*kx + c44*ky*ky + c33*kz*kz;
            a12= c11c66*kx*ky;
            a13= c13c44*kx*kz;
            a23= c13c44*ky*kz;

            Chr[0] = a11;
            Chr[4] = a22;
            Chr[8] = a33;
            Chr[1] = Chr[3] = a12;
            Chr[2] = Chr[6] = a13;
            Chr[5] = Chr[7] = a23;

            // LAPACK's ssyev routine (slow but accurate) 
            dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

            upx=Chr[6];
            upy=Chr[7];
            upz=Chr[8];

            if(upx*kx + upy*ky+ upz*kz < 0.) {
                upx=-Chr[6];
                upy=-Chr[7];
                upz=-Chr[8];
            }

            double usx= -kx*upz;
            double usy= -ky*upz;
            double usz= kx*upx + ky*upy;
 
            double rk=sqrt(usx*usx+usy*usy+usz*usz);
			if(rk==0.0)
               resz(a,b) = 0.0;
            else{
               usz /= rk;
			   double ss=sqrt(1.0-kz*kz);
               resz(a,b) = usz*ss;       /* scaling operator to tackle kiss singularity */
			}
              
         }// b loop
    }// a loop

    return 0;
}

static void map2d1d(float *d, DblNumMat mat, int m, int n)
{
   int i, j, k, kk;
   for(i=0; i < m; i++){
      kk = i*n;	  
      for (j=0; j < n; j++)
      {
        k = kk+j;
        d[k] = (float)mat(i,j);
	  }
   }

}

