/* 3-D two-components elastic wavefield extrapolation 
   using low-rank approximate PS solution on the base of 
   displacement wave equation in ORT media.

   Copyright (C) 2014 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
   modified by Peng Zou in 2015
     
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
#include<fftw3.h>
#include<omp.h>

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
#include "fwpvtielowrank.h"
#include "decomplowrank.h"
#include "vti2tti.h"
#include "eigen3x3.h"
}

static float *c11,*c12,*c13,*c14,*c15,*c16,
			      *c22,*c23,*c24,*c25,*c26,
				       *c33,*c34,*c35,*c36,
					        *c44,*c45,*c46,
						         *c55,*c56,
							          *c66;
static double dt1, dt2;

static std::valarray<double> rkx,rky,rkz,rk2;
static std::valarray<float> vp, vs, ep1, ep2, de1, de2, de3, ga1, ga2, th, ph;

/* dual-domain operators based on low-rank decomp. */
int sampleopx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopx3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopy1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopy2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopy3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopz3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

int sampleosx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosx3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosy1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosy2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosy3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosz3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

int sampleosz3(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

static void map2d1d(double *d, DblNumMat mat, int m, int n);

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);
   fftwf_init_threads();
   omp_set_num_threads(12);
   
   clock_t t1, t2, t3;
   float   timespent;

   t1=clock();

   int i,j,k;

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
   iRSF vp0, vs0("vs0"), epsi1("epsi1"),epsi2("epsi2"),del1("del1"),del2("del2"),
		del3("del3"), gam1("gam1"),gam2("gam2"),the("the"),phi("phi");

   /* Read/Write axes */
   int nxv,nyv,nzv;
   vp0.get("n1",nzv);
   vp0.get("n2",nxv);
   vp0.get("n3",nyv);

   float az, ax,ay;
   vp0.get("o1",az);
   vp0.get("o2",ax);
   vp0.get("o3",ay);

   float fx,fy,fz;
   fx=ax*1000.0;
   fy=ay*1000.0;
   fz=az*1000.0;

   float dx,dy,dz;
   vp0.get("d1",az);
   vp0.get("d2",ax);
   vp0.get("d3",ay);
   dz = az*1000.0;
   dx = ax*1000.0;
   dy = ay*1000.0;

   /* wave modeling space */
   int nx,ny,nz,nxz,nxyz;
   nx=nxv;
   ny=nyv;
   nz=nzv;
   nxz=nx*nz;
   nxyz=nx*ny*nz;

   sf_warning("nx=%d ny=%d nz=%d",nx,ny,nz);
   sf_warning("dx=%f dy=%f dz=%f",dx,dy,dz);

   sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("Warning: 2nd-order spectral need odd-based FFT");
   sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

   vp.resize(nxyz);
   vs.resize(nxyz);
   ep1.resize(nxyz);
   ep2.resize(nxyz);
   de1.resize(nxyz);
   de2.resize(nxyz);
   de3.resize(nxyz);
   ga1.resize(nxyz);
   ga2.resize(nxyz);
   th.resize(nxyz);
   ph.resize(nxyz);
 
   c11 = sf_floatalloc(nxyz);
   c12 = sf_floatalloc(nxyz);
   c13 = sf_floatalloc(nxyz);
   c14 = sf_floatalloc(nxyz);
   c15 = sf_floatalloc(nxyz);
   c16 = sf_floatalloc(nxyz);
   c22 = sf_floatalloc(nxyz);
   c23 = sf_floatalloc(nxyz);
   c24 = sf_floatalloc(nxyz);
   c25 = sf_floatalloc(nxyz);
   c26 = sf_floatalloc(nxyz);
   c33 = sf_floatalloc(nxyz);
   c34 = sf_floatalloc(nxyz);
   c35 = sf_floatalloc(nxyz);
   c36 = sf_floatalloc(nxyz);
   c44 = sf_floatalloc(nxyz);
   c45 = sf_floatalloc(nxyz);
   c46 = sf_floatalloc(nxyz);
   c55 = sf_floatalloc(nxyz);
   c56 = sf_floatalloc(nxyz);
   c66 = sf_floatalloc(nxyz);
   
   vp0>>vp;
   vs0>>vs;
   epsi1>>ep1;
   epsi2>>ep2;
   del1>>de1;
   del2>>de2;
   del3>>de3;
   gam1>>ga1;
   gam2>>ga2;
   the>>th;
   phi>>ph;

   for(i=0;i<nxyz;i++)
   {
       th[i] *= SF_PI/180.0;
       ph[i] *= SF_PI/180.0;
   }

   float *vp_1,*vs_1,*ep_1,*ep_2,*de_1,*de_2,*de_3,*ga_1,*ga_2,*th_1,*ph_1;
   vp_1 = sf_floatalloc(nxyz);
   vs_1 = sf_floatalloc(nxyz);
   ep_1 = sf_floatalloc(nxyz);
   ep_2 = sf_floatalloc(nxyz);
   de_1 = sf_floatalloc(nxyz);
   de_2 = sf_floatalloc(nxyz);
   de_3 = sf_floatalloc(nxyz);
   ga_1 = sf_floatalloc(nxyz);
   ga_2 = sf_floatalloc(nxyz);
   th_1 = sf_floatalloc(nxyz);
   ph_1 = sf_floatalloc(nxyz);

   for(int i=0;i<nxyz;i++)
   {
	   vp_1[i] = vp[i];
	   vs_1[i] = vs[i];
	   de_1[i] = de1[i];
	   de_2[i] = de2[i];
	   de_3[i] = de3[i];
	   ep_1[i] = ep1[i];
	   ep_2[i] = ep2[i];
	   ga_1[i] = ga1[i];
	   ga_2[i] = ga2[i];
	   th_1[i] = th[i];
	   ph_1[i] = ph[i];
   }
   Thomson2stiffness_ort(vp_1,vs_1,ep_1,ep_2,de_1,de_2,de_3,ga_1,ga_2,th_1,ph_1,c11,c12,c13,c14,c15,c16,
						c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,nx, ny, nz);

	free(vp_1);
	free(vs_1);
	free(ep_1);
	free(ep_2);
	free(de_1);
	free(de_2);
	free(de_3);
	free(ga_1);
	free(ga_2);
	free(th_1);
	free(ph_1);

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nkx*nky*nkz;

   float dkz,dkx,dky,kz0,kx0,ky0;

   dkx=2*SF_PI/dx/nx;
   dky=2*SF_PI/dy/ny;
   dkz=2*SF_PI/dz/nz;

   kx0=-SF_PI/dx;
   ky0=-SF_PI/dy;
   kz0=-SF_PI/dz;

   sf_warning("dkx=%f dky=%f dkz=%f",dkx,dky,dkz);

   rkx.resize(nk);
   rky.resize(nk);
   rkz.resize(nk);
   rk2.resize(nk);


   double kx, kz,ky,rk, k2;
   int ix,iy,iz;
   i = 0;
   for(iy=0; iy < nky; iy++)
   {
     ky = ky0+iy*dky;

     for(ix=0; ix < nkx; ix++)
     {
       kx = kx0+ix*dkx;

         for (iz=0; iz < nkz; iz++)
         {
            kz = kz0+iz*dkz;

            k2 = ky*ky+kx*kx+kz*kz;
            rk = sqrt(k2);

            rky[i] = ky/rk;
            rkx[i] = kx/rk;
            rkz[i] = kz/rk;
            rk2[i] = k2;
            i++;
         }
      }
   }

   vector<int> md(nxyz), nd(nk);
   for (k=0; k < nxyz; k++)  md[k] = k;
   for (k=0; k < nk; k++)  nd[k] = k;

   vector<int> lid, rid;
   DblNumMat mid, mat;

   int jump = 2;
   vector<int> ms, nss, js;
   ms.resize(3); ms[0] = nz; ms[1] = nx; ms[2] = ny;
   nss.resize(3); nss[0] = nz; nss[1] = nx; nss[2] = ny;
   js.resize(3); js[0] = jump;  js[1] = jump;  js[2] = jump;

   /********* qP-wave low rank decomposition of operator BxAx + BxyAxy + BxzAxz applying to ux for updating ux **********/
   int   m2opx1, n2opx1;
   double *ldataopx1, *fmidopx1, *rdataopx1;

   iC( ddlowrank(ms,nss,js,sampleopx1,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopx1,eps,npk,lid,rid,mid) );
   m2opx1=mid.m();
   n2opx1=mid.n();
   sf_warning("m2opx1=%d n2opx1=%d",m2opx1, n2opx1);

   fmidopx1  = (double*)malloc(sizeof(double)*m2opx1*n2opx1);
   ldataopx1 = (double*)malloc(sizeof(double)*nxyz*m2opx1);
   rdataopx1 = (double*)malloc(sizeof(double)*n2opx1*nk);

   map2d1d(fmidopx1, mid, m2opx1, n2opx1);

   iC ( sampleopx1(md,lid,mat) );
   map2d1d(ldataopx1, mat, nxyz, m2opx1);

   iC ( sampleopx1(rid,nd,mat) );
   map2d1d(rdataopx1, mat, n2opx1, nk);

   /********* qP-wave low rank decomposition of operator BxAxy + BxyAy + BxzAyz applying to uy for updating ux **********/
   int   m2opx2, n2opx2;
   double *ldataopx2, *fmidopx2, *rdataopx2;

   iC( ddlowrank(ms,nss,js,sampleopx2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopx2,eps,npk,lid,rid,mid) );
   m2opx2=mid.m();
   n2opx2=mid.n();
   sf_warning("m2opx2=%d n2opx2=%d",m2opx2, n2opx2);

   fmidopx2  = (double*)malloc(sizeof(double)*m2opx2*n2opx2);
   ldataopx2 = (double*)malloc(sizeof(double)*nxyz*m2opx2);
   rdataopx2 = (double*)malloc(sizeof(double)*n2opx2*nk);

   map2d1d(fmidopx2, mid, m2opx2, n2opx2);

   iC ( sampleopx2(md,lid,mat) );
   map2d1d(ldataopx2, mat, nxyz, m2opx2);

   iC ( sampleopx2(rid,nd,mat) );
   map2d1d(rdataopx2, mat, n2opx2, nk);

   /********* qP-wave low rank decomposition of operator BxAxz + BxyAyz + BxzAz applying to uz for updating ux **********/
   int   m2opx3, n2opx3;
   double *ldataopx3, *fmidopx3, *rdataopx3;

   iC( ddlowrank(ms,nss,js,sampleopx3,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopx3,eps,npk,lid,rid,mid) );
   m2opx3=mid.m();
   n2opx3=mid.n();
   sf_warning("m2opx3=%d n2opx3=%d",m2opx3, n2opx3);

   fmidopx3  = (double*)malloc(sizeof(double)*m2opx3*n2opx3);
   ldataopx3 = (double*)malloc(sizeof(double)*nxyz*m2opx3);
   rdataopx3 = (double*)malloc(sizeof(double)*n2opx3*nk);

   map2d1d(fmidopx3, mid, m2opx3, n2opx3);

   iC ( sampleopx3(md,lid,mat) );
   map2d1d(ldataopx3, mat, nxyz, m2opx3);

   iC ( sampleopx3(rid,nd,mat) );
   map2d1d(rdataopx3, mat, n2opx3, nk);

   /********* qP-wave low rank decomposition of operator BxyAx + ByAxy + ByzAxz applying to ux for updating uy **********/
   int   m2opy1, n2opy1;
   double *ldataopy1, *fmidopy1, *rdataopy1;

   iC( ddlowrank(ms,nss,js,sampleopy1,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopy1,eps,npk,lid,rid,mid) );
   m2opy1=mid.m();
   n2opy1=mid.n();
   sf_warning("m2opy1=%d n2opy1=%d",m2opy1, n2opy1);

   fmidopy1  = (double*)malloc(sizeof(double)*m2opy1*n2opy1);
   ldataopy1 = (double*)malloc(sizeof(double)*nxyz*m2opy1);
   rdataopy1 = (double*)malloc(sizeof(double)*n2opy1*nk);
   
   map2d1d(fmidopy1, mid, m2opy1, n2opy1);

   iC ( sampleopy1(md,lid,mat) );
   map2d1d(ldataopy1, mat, nxyz, m2opy1);

   iC ( sampleopy1(rid,nd,mat) );
   map2d1d(rdataopy1, mat, n2opy1, nk);
   
   /********* qP-wave low rank decomposition of operator BxyAxy + ByAy + ByzAyz applying to uy for updating uy **********/
   int   m2opy2, n2opy2;
   double *ldataopy2, *fmidopy2, *rdataopy2;
   
   iC( ddlowrank(ms,nss,js,sampleopy2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopy2,eps,npk,lid,rid,mid) );
   m2opy2=mid.m();
   n2opy2=mid.n();
   sf_warning("m2opy2=%d n2opy2=%d",m2opy2, n2opy2);

   fmidopy2  = (double*)malloc(sizeof(double)*m2opy2*n2opy2);
   ldataopy2 = (double*)malloc(sizeof(double)*nxyz*m2opy2);
   rdataopy2 = (double*)malloc(sizeof(double)*n2opy2*nk);
   
   map2d1d(fmidopy2, mid, m2opy2, n2opy2);

   iC ( sampleopy2(md,lid,mat) );
   map2d1d(ldataopy2, mat, nxyz, m2opy2);

   iC ( sampleopy2(rid,nd,mat) );
   map2d1d(rdataopy2, mat, n2opy2, nk);
   
   /********* qP-wave low rank decomposition of operator BxyAxz + ByAyz + ByzAz applying to uz for updating uy **********/
   int   m2opy3, n2opy3;
   double *ldataopy3, *fmidopy3, *rdataopy3;

   iC( ddlowrank(ms,nss,js,sampleopy3,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopy3,eps,npk,lid,rid,mid) );
   m2opy3=mid.m();
   n2opy3=mid.n();
   sf_warning("m2opy3=%d n2opy3=%d",m2opy3, n2opy3);

   fmidopy3  = (double*)malloc(sizeof(double)*m2opy3*n2opy3);
   ldataopy3 = (double*)malloc(sizeof(double)*nxyz*m2opy3);
   rdataopy3 = (double*)malloc(sizeof(double)*n2opy3*nk);

   map2d1d(fmidopy3, mid, m2opy3, n2opy3);

   iC ( sampleopy3(md,lid,mat) );
   map2d1d(ldataopy3, mat, nxyz, m2opy3);

   iC ( sampleopy3(rid,nd,mat) );
   map2d1d(rdataopy3, mat, n2opy3, nk);

   /********* qP-wave low rank decomposition of operator BxzAx + ByzAxy + BzAxz applying to ux for updating uz **********/
   int   m2opz1, n2opz1;
   double *ldataopz1, *fmidopz1, *rdataopz1;

   iC( ddlowrank(ms,nss,js,sampleopz1,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopz1,eps,npk,lid,rid,mid) );
   m2opz1=mid.m();
   n2opz1=mid.n();
   sf_warning("m2opz1=%d n2opz1=%d",m2opz1, n2opz1);

   fmidopz1  = (double*)malloc(sizeof(double)*m2opz1*n2opz1);
   ldataopz1 = (double*)malloc(sizeof(double)*nxyz*m2opz1);
   rdataopz1 = (double*)malloc(sizeof(double)*n2opz1*nk);

   map2d1d(fmidopz1, mid, m2opz1, n2opz1);

   iC ( sampleopz1(md,lid,mat) );
   map2d1d(ldataopz1, mat, nxyz, m2opz1);

   iC ( sampleopz1(rid,nd,mat) );
   map2d1d(rdataopz1, mat, n2opz1, nk);
   
   /********* qP-wave low rank decomposition of operator BxzAxy + ByzAy + BzAyz applying to uy for updating uz **********/
   int   m2opz2, n2opz2;
   double *ldataopz2, *fmidopz2, *rdataopz2;
   
   iC( ddlowrank(ms,nss,js,sampleopz2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopz2,eps,npk,lid,rid,mid) );
   m2opz2=mid.m();
   n2opz2=mid.n();
   sf_warning("m2opz2=%d n2opz2=%d",m2opz2, n2opz2);

   fmidopz2  = (double*)malloc(sizeof(double)*m2opz2*n2opz2);
   ldataopz2 = (double*)malloc(sizeof(double)*nxyz*m2opz2);
   rdataopz2 = (double*)malloc(sizeof(double)*n2opz2*nk);
   
   map2d1d(fmidopz2, mid, m2opz2, n2opz2);

   iC ( sampleopz2(md,lid,mat) );
   map2d1d(ldataopz2, mat, nxyz, m2opz2);

   iC ( sampleopz2(rid,nd,mat) );
   map2d1d(rdataopz2, mat, n2opz2, nk);
   
   /********* qP-wave low rank decomposition of operator BxzAxz + ByzAyz + BzAz applying to uz for updating uz **********/
   int   m2opz3, n2opz3;
   double *ldataopz3, *fmidopz3, *rdataopz3;
   
   iC( ddlowrank(ms,nss,js,sampleopz3,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleopz3,eps,npk,lid,rid,mid) );
   m2opz3=mid.m();
   n2opz3=mid.n();
   sf_warning("m2opz3=%d n2opz3=%d",m2opz3, n2opz3);

   fmidopz3  = (double*)malloc(sizeof(double)*m2opz3*n2opz3);
   ldataopz3 = (double*)malloc(sizeof(double)*nxyz*m2opz3);
   rdataopz3 = (double*)malloc(sizeof(double)*n2opz3*nk);

   map2d1d(fmidopz3, mid, m2opz3, n2opz3);

   iC ( sampleopz3(md,lid,mat) );
   map2d1d(ldataopz3, mat, nxyz, m2opz3);

   iC ( sampleopz3(rid,nd,mat) );
   map2d1d(rdataopz3, mat, n2opz3, nk);

   /********* qS-wave low rank decomposition of operator BxAx + BxyAxy + BxzAxz applying to ux for updating ux **********/
   int   m2osx1, n2osx1;
   double *ldataosx1, *fmidosx1, *rdataosx1;

   iC( ddlowrank(ms,nss,js,sampleosx1,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosx1,eps,npk,lid,rid,mid) );
   m2osx1=mid.m();
   n2osx1=mid.n();
   sf_warning("m2osx1=%d n2osx1=%d",m2osx1, n2osx1);

   fmidosx1  = (double*)malloc(sizeof(double)*m2osx1*n2osx1);
   ldataosx1 = (double*)malloc(sizeof(double)*nxyz*m2osx1);
   rdataosx1 = (double*)malloc(sizeof(double)*n2osx1*nk);

   map2d1d(fmidosx1, mid, m2osx1, n2osx1);

   iC ( sampleosx1(md,lid,mat) );
   map2d1d(ldataosx1, mat, nxyz, m2osx1);

   iC ( sampleosx1(rid,nd,mat) );
   map2d1d(rdataosx1, mat, n2osx1, nk);

   /********* qS-wave low rank decomposition of operator BxAxy + BxyAy + BxzAyz applying to uy for updating ux **********/
   int   m2osx2, n2osx2;
   double *ldataosx2, *fmidosx2, *rdataosx2;
   
   iC( ddlowrank(ms,nss,js,sampleosx2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosx2,eps,npk,lid,rid,mid) );
   m2osx2=mid.m();
   n2osx2=mid.n();
   sf_warning("m2osx2=%d n2osx2=%d",m2osx2, n2osx2);

   fmidosx2  = (double*)malloc(sizeof(double)*m2osx2*n2osx2);
   ldataosx2 = (double*)malloc(sizeof(double)*nxyz*m2osx2);
   rdataosx2 = (double*)malloc(sizeof(double)*n2osx2*nk);
   
   map2d1d(fmidosx2, mid, m2osx2, n2osx2);

   iC ( sampleosx2(md,lid,mat) );
   map2d1d(ldataosx2, mat, nxyz, m2osx2);

   iC ( sampleosx2(rid,nd,mat) );
   map2d1d(rdataosx2, mat, n2osx2, nk);
   
   /********* qS-wave low rank decomposition of operator BxAxz + BxyAyz + BxzAz applying to uz for updating ux **********/
   int   m2osx3, n2osx3;
   double *ldataosx3, *fmidosx3, *rdataosx3;
   
   iC( ddlowrank(ms,nss,js,sampleosx3,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosx3,eps,npk,lid,rid,mid) );
   m2osx3=mid.m();
   n2osx3=mid.n();
   sf_warning("m2osx3=%d n2osx3=%d",m2osx3, n2osx3);

   fmidosx3  = (double*)malloc(sizeof(double)*m2osx3*n2osx3);
   ldataosx3 = (double*)malloc(sizeof(double)*nxyz*m2osx3);
   rdataosx3 = (double*)malloc(sizeof(double)*n2osx3*nk);
   
   map2d1d(fmidosx3, mid, m2osx3, n2osx3);

   iC ( sampleosx3(md,lid,mat) );
   map2d1d(ldataosx3, mat, nxyz, m2osx3);

   iC ( sampleosx3(rid,nd,mat) );
   map2d1d(rdataosx3, mat, n2osx3, nk);

   /********* qS-wave low rank decomposition of operator BxyAx + ByAxy + ByzAxz applying to ux for updating uy **********/
   int   m2osy1, n2osy1;
   double *ldataosy1, *fmidosy1, *rdataosy1;

   iC( ddlowrank(ms,nss,js,sampleosy1,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosy1,eps,npk,lid,rid,mid) );
   m2osy1=mid.m();
   n2osy1=mid.n();
   sf_warning("m2osy1=%d n2osy1=%d",m2osy1, n2osy1);

   fmidosy1  = (double*)malloc(sizeof(double)*m2osy1*n2osy1);
   ldataosy1 = (double*)malloc(sizeof(double)*nxyz*m2osy1);
   rdataosy1 = (double*)malloc(sizeof(double)*n2osy1*nk);

   map2d1d(fmidosy1, mid, m2osy1, n2osy1);

   iC ( sampleosy1(md,lid,mat) );
   map2d1d(ldataosy1, mat, nxyz, m2osy1);

   iC ( sampleosy1(rid,nd,mat) );
   map2d1d(rdataosy1, mat, n2osy1, nk);

   /********* qS-wave low rank decomposition of operator BxyAxy + ByAy + ByzAyz applying to uy for updating uy **********/
   int   m2osy2, n2osy2;
   double *ldataosy2, *fmidosy2, *rdataosy2;

   iC( ddlowrank(ms,nss,js,sampleosy2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosy2,eps,npk,lid,rid,mid) );
   m2osy2=mid.m();
   n2osy2=mid.n();
   sf_warning("m2osy2=%d n2osy2=%d",m2osy2, n2osy2);

   fmidosy2  = (double*)malloc(sizeof(double)*m2osy2*n2osy2);
   ldataosy2 = (double*)malloc(sizeof(double)*nxyz*m2osy2);
   rdataosy2 = (double*)malloc(sizeof(double)*n2osy2*nk);
   
   map2d1d(fmidosy2, mid, m2osy2, n2osy2);

   iC ( sampleosy2(md,lid,mat) );
   map2d1d(ldataosy2, mat, nxyz, m2osy2);

   iC ( sampleosy2(rid,nd,mat) );
   map2d1d(rdataosy2, mat, n2osy2, nk);
   
   /********* qS-wave low rank decomposition of operator BxyAxz + ByAyz + ByzAz applying to uz for updating uy **********/
   int   m2osy3, n2osy3;
   double *ldataosy3, *fmidosy3, *rdataosy3;
   
   iC( ddlowrank(ms,nss,js,sampleosy3,eps,npk,lid,rid,mid) );
  // iC( ddlowrank(nxyz,nk,sampleosy3,eps,npk,lid,rid,mid) );
   m2osy3=mid.m();
   n2osy3=mid.n();
   sf_warning("m2osy3=%d n2osy3=%d",m2osy3, n2osy3);

   fmidosy3  = (double*)malloc(sizeof(double)*m2osy3*n2osy3);
   ldataosy3 = (double*)malloc(sizeof(double)*nxyz*m2osy3);
   rdataosy3 = (double*)malloc(sizeof(double)*n2osy3*nk);
   
   map2d1d(fmidosy3, mid, m2osy3, n2osy3);

   iC ( sampleosy3(md,lid,mat) );
   map2d1d(ldataosy3, mat, nxyz, m2osy3);

   iC ( sampleosy3(rid,nd,mat) );
   map2d1d(rdataosy3, mat, n2osy3, nk);
   
   /********* qS-wave low rank decomposition of operator BxzAx + ByzAxy + BzAxz applying to ux for updating uz **********/
   int   m2osz1, n2osz1;
   double *ldataosz1, *fmidosz1, *rdataosz1;

   iC( ddlowrank(ms,nss,js,sampleosz1,eps,npk,lid,rid,mid) );
  // iC( ddlowrank(nxyz,nk,sampleosz1,eps,npk,lid,rid,mid) );
   m2osz1=mid.m();
   n2osz1=mid.n();
   sf_warning("m2osz1=%d n2osz1=%d",m2osz1, n2osz1);

   fmidosz1  = (double*)malloc(sizeof(double)*m2osz1*n2osz1);
   ldataosz1 = (double*)malloc(sizeof(double)*nxyz*m2osz1);
   rdataosz1 = (double*)malloc(sizeof(double)*n2osz1*nk);

   map2d1d(fmidosz1, mid, m2osz1, n2osz1);

   iC ( sampleosz1(md,lid,mat) );
   map2d1d(ldataosz1, mat, nxyz, m2osz1);

   iC ( sampleosz1(rid,nd,mat) );
   map2d1d(rdataosz1, mat, n2osz1, nk);

   /********* qS-wave low rank decomposition of operator BxzAxy + ByzAy + BzAyz applying to uy for updating uz **********/
   int   m2osz2, n2osz2;
   double *ldataosz2, *fmidosz2, *rdataosz2;

   iC( ddlowrank(ms,nss,js,sampleosz2,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosz2,eps,npk,lid,rid,mid) );
   m2osz2=mid.m();
   n2osz2=mid.n();
   sf_warning("m2osz2=%d n2osz2=%d",m2osz2, n2osz2);

   fmidosz2  = (double*)malloc(sizeof(double)*m2osz2*n2osz2);
   ldataosz2 = (double*)malloc(sizeof(double)*nxyz*m2osz2);
   rdataosz2 = (double*)malloc(sizeof(double)*n2osz2*nk);

   map2d1d(fmidosz2, mid, m2osz2, n2osz2);

   iC ( sampleosz2(md,lid,mat) );
   map2d1d(ldataosz2, mat, nxyz, m2osz2);

   iC ( sampleosz2(rid,nd,mat) );
   map2d1d(rdataosz2, mat, n2osz2, nk);
   
   /********* qS-wave low rank decomposition of operator BxzAxz + ByzAyz + BzAz applying to uz for updating uz **********/
   int   m2osz3, n2osz3;
   double *ldataosz3, *fmidosz3, *rdataosz3;
   
   iC( ddlowrank(ms,nss,js,sampleosz3,eps,npk,lid,rid,mid) );
   //iC( ddlowrank(nxyz,nk,sampleosz3,eps,npk,lid,rid,mid) );
   m2osz3=mid.m();
   n2osz3=mid.n();
   sf_warning("m2osz3=%d n2osz3=%d",m2osz3, n2osz3);

   fmidosz3  = (double*)malloc(sizeof(double)*m2osz3*n2osz3);
   ldataosz3 = (double*)malloc(sizeof(double)*nxyz*m2osz3);
   rdataosz3 = (double*)malloc(sizeof(double)*n2osz3*nk);
   
   map2d1d(fmidosz3, mid, m2osz3, n2osz3);

   iC ( sampleosz3(md,lid,mat) );
   map2d1d(ldataosz3, mat, nxyz, m2osz3);

   iC ( sampleosz3(rid,nd,mat) );
   map2d1d(rdataosz3, mat, n2osz3, nk);

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

   sf_warning("fx=%f fy=%f fz=%f ",fx,fy,fz);
   sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);

   /* source definition */
   int ixs, izs,iys;
   ixs=nxv/2;
   izs=nxv/2;
   iys=nyv/2;
   sf_warning("ixs=%d iys=%d izs=%d ", ixs,iys,izs);

   /* setup I/O files */
   oRSF Elasticx("out");
   oRSF Elasticy("Elasticy");
   oRSF Elasticz("Elasticz");
   oRSF ElasticPx("ElasticPx");
   oRSF ElasticPy("ElasticPy");
   oRSF ElasticPz("ElasticPz");
   oRSF ElasticSx("ElasticSx");
   oRSF ElasticSy("ElasticSy");
   oRSF ElasticSz("ElasticSz");

   Elasticx.put("n1",nz);
   Elasticx.put("n2",nx);
   Elasticx.put("n3",ny);
   Elasticx.put("d1",dz/1000);
   Elasticx.put("d2",dx/1000);
   Elasticx.put("d3",dy/1000);
   Elasticx.put("o1",fz/1000);
   Elasticx.put("o2",fx/1000);
   Elasticx.put("o3",fy/1000);

   Elasticy.put("n1",nz);
   Elasticy.put("n2",nx);
   Elasticy.put("n3",ny);
   Elasticy.put("d1",dz/1000);
   Elasticy.put("d2",dx/1000);
   Elasticy.put("o1",fz/1000);
   Elasticy.put("o2",fx/1000);
   
   Elasticz.put("n1",nkz);
   Elasticz.put("n2",nkx);
   Elasticz.put("d1",dz/1000);
   Elasticz.put("d2",dx/1000);
   Elasticz.put("d3",dy/1000);
   Elasticz.put("o1",fz/1000);
   Elasticz.put("o2",fx/1000);
   Elasticz.put("o3",fy/1000);
   
   ElasticPx.put("n1",nz);
   ElasticPx.put("n2",nx);
   ElasticPx.put("n3",ny);
   ElasticPx.put("d1",dz/1000);
   ElasticPx.put("d2",dx/1000);
   ElasticPx.put("d3",dy/1000);
   ElasticPx.put("o1",fz/1000);
   ElasticPx.put("o2",fx/1000);
   ElasticPx.put("o3",fy/1000);
   
   ElasticPy.put("n1",nz);
   ElasticPy.put("n2",nx);
   ElasticPy.put("n3",ny);
   ElasticPy.put("d1",dz/1000);
   ElasticPy.put("d2",dx/1000);
   ElasticPy.put("d3",dy/1000);
   ElasticPy.put("o1",fz/1000);
   ElasticPy.put("o2",fx/1000);
   ElasticPy.put("o3",fy/1000);

   ElasticPz.put("n1",nz);
   ElasticPz.put("n2",nx);
   ElasticPz.put("n3",ny);
   ElasticPz.put("d1",dz/1000);
   ElasticPz.put("d2",dx/1000);
   ElasticPz.put("d3",dy/1000);
   ElasticPz.put("o1",fz/1000);
   ElasticPz.put("o2",fx/1000);
   ElasticPz.put("o3",fy/1000);

   ElasticSx.put("n1",nz);
   ElasticSx.put("n2",nx);
   ElasticSx.put("n3",ny);
   ElasticSx.put("d1",dz/1000);
   ElasticSx.put("d2",dx/1000);
   ElasticSx.put("d3",dy/1000);
   ElasticSx.put("o1",fz/1000);
   ElasticSx.put("o2",fx/1000);
   ElasticSx.put("o3",fy/1000);

   ElasticSy.put("n1",nz);
   ElasticSy.put("n2",nx);
   ElasticSy.put("n3",ny);
   ElasticSy.put("d1",dz/1000);
   ElasticSy.put("d2",dx/1000);
   ElasticSy.put("d3",dy/1000);
   ElasticSy.put("o1",fz/1000);
   ElasticSy.put("o2",fx/1000);
   ElasticSy.put("o3",fy/1000);

   ElasticSz.put("n1",nz);
   ElasticSz.put("n2",nx);
   ElasticSz.put("n3",ny);
   ElasticSz.put("d1",dz/1000);
   ElasticSz.put("d2",dx/1000);
   ElasticSz.put("d3",dy/1000);
   ElasticSz.put("o1",fz/1000);
   ElasticSz.put("o2",fx/1000);
   ElasticSz.put("o3",fy/1000);

   /********************* wavefield extrapolation *************************/
   double *px1=(double*)malloc(sizeof(double)*nxyz);
   double *px2=(double*)malloc(sizeof(double)*nxyz);
   double *px3=(double*)malloc(sizeof(double)*nxyz);
   double *py1=(double*)malloc(sizeof(double)*nxyz);
   double *py2=(double*)malloc(sizeof(double)*nxyz);
   double *py3=(double*)malloc(sizeof(double)*nxyz);
   double *pz1=(double*)malloc(sizeof(double)*nxyz);
   double *pz2=(double*)malloc(sizeof(double)*nxyz);
   double *pz3=(double*)malloc(sizeof(double)*nxyz);
   
   double *sx1=(double*)malloc(sizeof(double)*nxyz);
   double *sx2=(double*)malloc(sizeof(double)*nxyz);
   double *sx3=(double*)malloc(sizeof(double)*nxyz);
   double *sy1=(double*)malloc(sizeof(double)*nxyz);
   double *sy2=(double*)malloc(sizeof(double)*nxyz);
   double *sy3=(double*)malloc(sizeof(double)*nxyz);
   double *sz1=(double*)malloc(sizeof(double)*nxyz);
   double *sz2=(double*)malloc(sizeof(double)*nxyz);
   double *sz3=(double*)malloc(sizeof(double)*nxyz);
   
   double *ux=(double*)malloc(sizeof(double)*nxyz);
   double *uy=(double*)malloc(sizeof(double)*nxyz);
   double *uz=(double*)malloc(sizeof(double)*nxyz);

   double *pp=(double*)malloc(sizeof(double)*nxyz);
   double *ppp=(double*)malloc(sizeof(double)*nxyz);
   
   zero1double(px1, nxyz);
   zero1double(px2, nxyz);
   zero1double(px3, nxyz);
   zero1double(py1, nxyz);
   zero1double(py2, nxyz);
   zero1double(py3, nxyz);
   zero1double(pz1, nxyz);
   zero1double(pz2, nxyz);
   zero1double(pz3, nxyz);

   zero1double(sx1, nxyz);
   zero1double(sx2, nxyz);
   zero1double(sx3, nxyz);
   zero1double(sy1, nxyz);
   zero1double(sy2, nxyz);
   zero1double(sy3, nxyz);
   zero1double(sz1, nxyz);
   zero1double(sz2, nxyz);
   zero1double(sz3, nxyz);

   zero1double(ux, nxyz);
   zero1double(uy, nxyz);
   zero1double(uz, nxyz);

   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);
   int *ijky = sf_intalloc(nky);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

   std::valarray<float> x(nxyz);

	int nbd = 30;
	double alpha = 0.001;
    double *decay = (double*)malloc(sizeof(double)*nxyz);
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            for (iz=0; iz < nz; iz++) {
                i = iz+nz *(ix+nx *iy);
                decay[i]=1.0;
                if(iz<nbd)
                    decay[i] *= exp(-pow(alpha*(nbd-iz)*dz,2));
                else if(iz>(nz-1-nbd))
                    decay[i] *= exp(-pow(alpha*(iz-nz+nbd+1)*dz,2));
                if(ix<nbd)
                    decay[i] *= exp(-pow(alpha*(nbd-ix)*dx,2));
                else if(ix>(nx-1-nbd))
                    decay[i] *= exp(-pow(alpha*(ix-nx+nbd+1)*dx,2));
                if(iy<nbd)
                    decay[i] *= exp(-pow(alpha*(nbd-iy)*dy,2));
                else if(iy>(ny-1-nbd))
                    decay[i] *= exp(-pow(alpha*(iy-ny+nbd+1)*dy,2));
            }
        }
    }

   int iii;
   for(int it=0;it<ns;it++)
   {
        float t=it*dt;

        if(it%10==0)
                sf_warning("Elastic: it= %d  t=%f(s)",it,t);
 
         // 3D exploding force source (e.g., Wu's PhD)
         for(k=-1;k<=1;k++)
         for(i=-1;i<=1;i++)
         for(j=-1;j<=1;j++)
         {
             if(fabs(k)+fabs(i)+fabs(j)==3)
             {
                 iii=(iys+k)*nxz+(ixs+i)*nz+(izs+j);
                 uy[iii]+=k*Ricker(t, f0, t0, A);
                 ux[iii]+=i*Ricker(t, f0, t0, A);
                 uz[iii]+=j*Ricker(t, f0, t0, A);
             }
        }
        // 3D 45-degree force source 
        //p2[isy][isx][isz]+=Ricker(t, f0, t0, A);
        //q2[isy][isx][isz]+=Ricker(t, f0, t0, A);
        //r2[isy][isx][isz]+=Ricker(t, f0, t0, A);

        if(it%10==0) sf_warning("ux=%f uy=%f uz=%f ",ux[iii],uy[iii],uz[iii]);

       /* extrapolation of Upx-componet */
        fwpvti3delowrank_double(ldataopx1,rdataopx1,fmidopx1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opx1,n2opx1);
        for(i=0;i<nxyz;i++) px3[i] = pp[i];
        fwpvti3delowrank_double(ldataopx2,rdataopx2,fmidopx2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opx2,n2opx2);
        for(i=0;i<nxyz;i++) px3[i] += pp[i];
        fwpvti3delowrank_double(ldataopx3,rdataopx3,fmidopx3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opx3,n2opx3);
        for(i=0;i<nxyz;i++) px3[i] = decay[i]*(px3[i] + pp[i]) - pow(decay[i],2)*px1[i];

        /* extrapolation of Upy-componet */
        fwpvti3delowrank_double(ldataopy1,rdataopy1,fmidopy1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opy1,n2opy1);
        for(i=0;i<nxyz;i++) py3[i] = pp[i];
        fwpvti3delowrank_double(ldataopy2,rdataopy2,fmidopy2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opy2,n2opy2);
        for(i=0;i<nxyz;i++) py3[i] += pp[i];
        fwpvti3delowrank_double(ldataopy3,rdataopy3,fmidopy3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opy3,n2opy3);
        for(i=0;i<nxyz;i++) py3[i] = decay[i]*(py3[i] + pp[i]) - pow(decay[i],2)*py1[i];

        /* extrapolation of Upz-componet */
        fwpvti3delowrank_double(ldataopz1,rdataopz1,fmidopz1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opz1,n2opz1);
        for(i=0;i<nxyz;i++) pz3[i] = pp[i];
        fwpvti3delowrank_double(ldataopz2,rdataopz2,fmidopz2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opz2,n2opz2);
        for(i=0;i<nxyz;i++) pz3[i] += pp[i];
        fwpvti3delowrank_double(ldataopz3,rdataopz3,fmidopz3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2opz3,n2opz3);
        for(i=0;i<nxyz;i++) pz3[i] = decay[i]*(pz3[i] + pp[i]) - pow(decay[i],2)*pz1[i];

        /* extrapolation of Usx-componet */
        fwpvti3delowrank_double(ldataosx1,rdataosx1,fmidosx1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osx1,n2osx1);
        for(i=0;i<nxyz;i++) sx3[i] = pp[i];
        fwpvti3delowrank_double(ldataosx2,rdataosx2,fmidosx2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osx2,n2osx2);
        for(i=0;i<nxyz;i++) sx3[i] += pp[i];
        fwpvti3delowrank_double(ldataosx3,rdataosx3,fmidosx3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osx3,n2osx3);
        for(i=0;i<nxyz;i++) sx3[i] = decay[i]*(sx3[i] + pp[i]) - pow(decay[i],2)*sx1[i];

        /* extrapolation of Usy-componet */
        fwpvti3delowrank_double(ldataosy1,rdataosy1,fmidosy1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osy1,n2osy1);
        for(i=0;i<nxyz;i++) sy3[i] = pp[i];
        fwpvti3delowrank_double(ldataosy2,rdataosy2,fmidosy2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osy2,n2osy2);
        for(i=0;i<nxyz;i++) sy3[i] += pp[i];
        fwpvti3delowrank_double(ldataosy3,rdataosy3,fmidosy3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osy3,n2osy3);
        for(i=0;i<nxyz;i++) sy3[i] = decay[i]*(sy3[i] + pp[i]) - pow(decay[i],2)*sy1[i];

        /* extrapolation of Usz-componet */
        fwpvti3delowrank_double(ldataosz1,rdataosz1,fmidosz1,pp,ux,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osz1,n2osz1);
        for(i=0;i<nxyz;i++) sz3[i] = pp[i];
        fwpvti3delowrank_double(ldataosz2,rdataosz2,fmidosz2,pp,uy,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osz2,n2osz2);
        for(i=0;i<nxyz;i++) sz3[i] += pp[i];
        fwpvti3delowrank_double(ldataosz3,rdataosz3,fmidosz3,pp,uz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2osz3,n2osz3);
        for(i=0;i<nxyz;i++) sz3[i] = decay[i]*(sz3[i] + pp[i]) - pow(decay[i],2)*sz1[i];

        /******* update the wavefield ********/
        for(i=0;i<nxyz;i++){
                px1[i]=px2[i];
                px2[i]=px3[i];
                py1[i]=py2[i];
                py2[i]=py3[i];
                pz1[i]=pz2[i];
                pz2[i]=pz3[i];
                sx1[i]=sx2[i];
                sx2[i]=sx3[i];
                sy1[i]=sy2[i];
                sy2[i]=sy3[i];
                sz1[i]=sz2[i];
                sz2[i]=sz3[i];
                ux[i] = px3[i]+sx3[i];
                uy[i] = py3[i]+sy3[i];
                uz[i] = pz3[i]+sz3[i];
            }

        /******* output wavefields: components******/
        if(it==ns-1)
        {
              for(i=0;i<nxyz;i++) x[i]=ux[i];
              Elasticx<<x;
              for(i=0;i<nxyz;i++) x[i]=uy[i];
              Elasticy<<x;
              for(i=0;i<nxyz;i++) x[i]=uz[i];
              Elasticz<<x;
              for(i=0;i<nxyz;i++) x[i]=px3[i];
              ElasticPx<<x;
              for(i=0;i<nxyz;i++) x[i]=py3[i];
              ElasticPy<<x;
              for(i=0;i<nxyz;i++) x[i]=pz3[i];
              ElasticPz<<x;
              for(i=0;i<nxyz;i++) x[i]=sx3[i];
              ElasticSx<<x;
              for(i=0;i<nxyz;i++) x[i]=sy3[i];
              ElasticSy<<x;
              for(i=0;i<nxyz;i++) x[i]=sz3[i];
              ElasticSz<<x;
        }
   } //* it loop */

   t3=clock();
   timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
   sf_warning("CPU time for wavefield extrapolation.: %f(second)",timespent);

   timespent=(float)(t3-t1)/(ns*CLOCKS_PER_SEC);
   sf_warning("CPU time for every time extrapolation (including low-rank decom.): %f(second)",timespent);

   free(ldataopx1);
   free(fmidopx1);
   free(rdataopx1);

   free(ldataopx2);
   free(fmidopx2);
   free(rdataopx2);

   free(ldataopx3);
   free(fmidopx3);
   free(rdataopx3);

   free(ldataopy1);
   free(fmidopy1);
   free(rdataopy1);

   free(ldataopy2);
   free(fmidopy2);
   free(rdataopy2);

   free(ldataopy3);
   free(fmidopy3);
   free(rdataopy3);

   free(ldataopz1);
   free(fmidopz1);
   free(rdataopz1);

   free(ldataopz2);
   free(fmidopz2);
   free(rdataopz2);

   free(ldataopz3);
   free(fmidopz3);
   free(rdataopz3);

   free(ldataosx1);
   free(fmidosx1);
   free(rdataosx1);

   free(ldataosx2);
   free(fmidosx2);
   free(rdataosx2);

   free(ldataosx3);
   free(fmidosx3);
   free(rdataosx3);

   free(ldataosy1);
   free(fmidosy1);
   free(rdataosy1);

   free(ldataosy2);
   free(fmidosy2);
   free(rdataosy2);

   free(ldataosy3);
   free(fmidosy3);
   free(rdataosy3);

   free(ldataosz1);
   free(fmidosz1);
   free(rdataosz1);

   free(ldataosz2);
   free(fmidosz2);
   free(rdataosz2);

   free(ldataosz3);
   free(fmidosz3);
   free(rdataosz3);

   free(px1);
   free(px2);
   free(px3);
   free(py1);
   free(py2);
   free(py3);
   free(pz1);
   free(pz2);
   free(pz3);
   
   free(sx1);
   free(sx2);
   free(sx3);
   free(sy1);
   free(sy2);
   free(sy3);
   free(sz1);
   free(sz2);
   free(sz3);

   free(pp);
   free(ppp);

   free(ux);
   free(uy);
   free(uz);

   free(ijkx);
   free(ijky);
   free(ijkz);

   exit(0);
}

double A[3][3],Q[3][3],w[3]; //Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
int info;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xx*w_xx+d_xy*w_xy+d_xz*w_xz applying to ux **********//
int sampleopx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = u1*u1*axx + u1*u2*axy + u1*u3*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xx*w_xy+d_xy*w_yy+d_xz*w_yz applying to uy **********/
int sampleopx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = u1*u1*axy + u1*u2*ayy + u1*u3*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xx*w_xz+d_xy*w_yz+dxz*w_zz applying to uz **********/
int sampleopx3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = u1*u1*axz + u1*u2*ayz + u1*u3*azz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xy*w_xx+d_yy*w_xy+d_yz*wxz applying to ux **********//
int sampleopy1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = u1*u2*axx + u2*u2*axy + u2*u3*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xy*w_xy+d_yy*w_yy+d_yz*w_yz applying to uy **********/
int sampleopy2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = u2*u1*axy + u2*u2*ayy + u2*u3*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xy*w_xz+d_yy*w_yz+d_yz*w_zz applying to uz **********/
int sampleopy3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = u2*u1*axz + u2*u2*ayz + u2*u3*azz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xz*w_xx+d_yz*w_xy+d_yz*w_xz  applying to ux **********//
int sampleopz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = u3*u1*axx + u3*u2*axy + u3*u3*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xz*w_xy+d_yz*w_yy+d_zz*w_yz applying to uy **********/
int sampleopz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = u3*u1*axy + u3*u2*ayy + u3*u3*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qP-wave low rank decomposition of operator d_xz*w_xz+d_yz*w_yz+d_zz*w_zz applying to uz **********/
int sampleopz3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = u3*u1*axz + u3*u2*ayz + u3*u3*azz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xx*w_xx+d_xy*w_xy+d_xz*w_xz applying to ux **********//
int sampleosx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = (u2*u2+u3*u3)*axx - u1*u2*axy - u1*u3*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xx*w_xy+d_xy*w_yy+d_xz*w_yz applying to uy **********/
int sampleosx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = (u2*u2+u3*u3)*axy - u1*u2*ayy - u1*u3*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xx*w_xz+d_xy*w_yz+dxz*w_zz applying to uz **********/
int sampleosx3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = (u2*u2+u3*u3)*axz - u1*u2*ayz - u1*u3*azz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xy*w_xx+d_yy*w_xy+d_yz*wxz applying to ux **********//
int sampleosy1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = -u1*u2*axx + (u1*u1+u3*u3)*axy - u2*u3*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xy*w_xy+d_yy*w_yy+d_yz*w_yz applying to uy **********/
int sampleosy2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = -u2*u1*axy + (u1*u1+u3*u3)*ayy - u2*u3*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xy*w_xz+d_yy*w_yz+d_yz*w_zz applying to uz **********/
int sampleosy3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = -u2*u1*axz + (u1*u1+u3*u3)*ayz - u2*u3*azz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xz*w_xx+d_yz*w_xy+d_yz*w_xz  applying to ux **********//
int sampleosz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double axx = 2.0 - dt2*a11; 
			double axy = -dt2*a12;
			double axz = -dt2*a13;

            resx(a,b) = -u3*u1*axx - u3*u2*axy + (u1*u1+u2*u2)*axz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xz*w_xy+d_yz*w_yy+d_zz*w_yz applying to uy **********/
int sampleosz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double ayy = 2.0 - dt2*a22; 
			double axy = -dt2*a12;
			double ayz = -dt2*a23;

            resx(a,b) = -u3*u1*axy - u3*u2*ayy + (u1*u1+u2*u2)*ayz; 

         }// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//********* qS-wave low rank decomposition of operator d_xz*w_xz+d_yz*w_yz+d_zz*w_zz applying to uz **********/
int sampleosz3(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

	double a11, a12, a22, a33, a13, a23;
	double u1, u2, u3;
	double lam1,lam2,lam3,sinclam1,sinclam2,sinclam3;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];

        for(int b=0; b<nc; b++)
        {
            double kx = rkx[cs[b]];
            double ky = rky[cs[b]];
            double kz = rkz[cs[b]];
            double k2 = rk2[cs[b]];
            if(kx==0.0&&ky==0.0&&kz==0.0)
            {
               resx(a,b) = 0.0;
               continue;
            }

           double kx2 = kx*kx*k2;
           double ky2 = ky*ky*k2;
           double kz2 = kz*kz*k2;
		   double kxky = kx*ky*k2;
		   double kxkz = kx*kz*k2;
		   double kykz = ky*kz*k2;

		   a11 = c11[i]*kx2 + c66[i]*ky2 + c55[i]*kz2 + 2.0*(c56[i]*kykz + c15[i]*kxkz + c16[i]*kxky);
		   a22 = c66[i]*kx2 + c22[i]*ky2 + c44[i]*kz2 + 2.0*(c24[i]*kykz + c46[i]*kxkz + c26[i]*kxky);
		   a33 = c55[i]*kx2 + c44[i]*ky2 + c33[i]*kz2 + 2.0*(c34[i]*kykz + c35[i]*kxkz + c45[i]*kxky);
		   a12 = c16[i]*kx2 + c26[i]*ky2 + c45[i]*kz2 + (c46[i]+c25[i])*kykz + (c14[i]+c56[i])*kxkz + (c12[i]+c66[i])*kxky;
		   a13 = c15[i]*kx2 + c46[i]*ky2 + c35[i]*kz2 + (c45[i]+c36[i])*kykz + (c13[i]+c55[i])*kxkz + (c14[i]+c56[i])*kxky;
		   a23 = c56[i]*kx2 + c24[i]*ky2 + c34[i]*kz2 + (c44[i]+c23[i])*kykz + (c36[i]+c45[i])*kxkz + (c25[i]+c46[i])*kxky;

		   A[0][0] = a11;
		   A[0][1] = a12;
		   A[0][2] = a13;
		   A[1][0] = A[0][1];
		   A[1][1] = a22;
		   A[1][2] = a23;
		   A[2][0] = A[0][2];
		   A[2][1] = A[1][2];
		   A[2][2] = a33;

		   info = dsyevd3(A,Q,w);
		   if(info == -1)
		   {
			   sf_warning("Error in Calculation the eigenvalues and normalized eigenvectors");
			   exit(0);
		   }
		   u1 = Q[0][0];
		   u2 = Q[1][0];
		   u3 = Q[2][0];

		   if(u1*kx + u2*ky+ u3*kz < 0.) {
			   u1 = -u1;
			   u2 = -u2;
			   u3 = -u3;
		   }

		   lam1 = sqrt(w[0])*0.5*dt1;
		   lam2 = sqrt(w[1])*0.5*dt1;
		   lam3 = sqrt(w[2])*0.5*dt1;
		   sinclam1 = sin(lam1)*sin(lam1)/lam1/lam1;
		   sinclam2 = sin(lam2)*sin(lam2)/lam2/lam2;
		   sinclam3 = sin(lam3)*sin(lam3)/lam3/lam3;

		   a11 = Q[0][0]*Q[0][0]*w[0]*sinclam1 + Q[0][1]*Q[0][1]*w[1]*sinclam2 + Q[0][2]*Q[0][2]*w[2]*sinclam3;
		   a12 = Q[0][0]*Q[1][0]*w[0]*sinclam1 + Q[0][1]*Q[1][1]*w[1]*sinclam2 + Q[0][2]*Q[1][2]*w[2]*sinclam3;
		   a13 = Q[0][0]*Q[2][0]*w[0]*sinclam1 + Q[0][1]*Q[2][1]*w[1]*sinclam2 + Q[0][2]*Q[2][2]*w[2]*sinclam3;
		   a22 = Q[1][0]*Q[1][0]*w[0]*sinclam1 + Q[1][1]*Q[1][1]*w[1]*sinclam2 + Q[1][2]*Q[1][2]*w[2]*sinclam3;
		   a23 = Q[1][0]*Q[2][0]*w[0]*sinclam1 + Q[1][1]*Q[2][1]*w[1]*sinclam2 + Q[1][2]*Q[2][2]*w[2]*sinclam3;
		   a33 = Q[2][0]*Q[2][0]*w[0]*sinclam1 + Q[2][1]*Q[2][1]*w[1]*sinclam2 + Q[2][2]*Q[2][2]*w[2]*sinclam3;

           // wavefield extrapolator
			double azz = 2.0 - dt2*a33; 
			double axz = -dt2*a13;
			double ayz = -dt2*a23;

            resx(a,b) = -u3*u1*axz - u3*u2*ayz + (u1*u1+u2*u2)*azz; 

         }// b loop
    }// a loop

    return 0;
}

static void map2d1d(double *d, DblNumMat mat, int m, int n)
{
   int i, j, k;
   k=0;
   for (i=0; i < m; i++)
   for (j=0; j < n; j++)
   {
        d[k] = (double)mat(i,j);
        k++;
   }

}
