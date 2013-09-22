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

using namespace std;

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "kykxkztaper.h"
#include "sepdivcurl.h"
#include "puthead.h"
}

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4, t5, t44;
   float   timespent;

   t1=clock();

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

   sf_warning("fx=%f fy=%f fz=%f dx=%f dy=%f dz=%f",fx,fy,fz,dx,dy,dz);
   sf_warning("nx=%d ny=%d nz=%d ", nx,ny,nz);

   /* wave modeling space */
   int nxyz=nx*ny*nz;

   /* Fourier spectra demension */
   int nkz,nkx,nky,nk;
   nkx=nx;
   nky=ny;
   nkz=nz;
   nk = nky*nkx*nkz;

   sf_warning("nxyz=%d nk=%d",nxyz,nk);

   double dkz,dkx,dky,kz0,kx0,ky0;

   dkx=2*PI/dx/nx;
   dky=2*PI/dy/ny;
   dkz=2*PI/dz/nz;

   kx0=-PI/dx;
   ky0=-PI/dy;
   kz0=-PI/dz;

   double *rkx, *rky, *rkz;

   rkx=(double*)calloc(sizeof(double),nxyz);
   rky=(double*)calloc(sizeof(double),nxyz);
   rkz=(double*)calloc(sizeof(double),nxyz);

   double kx, ky, kz, k2, rk;
   int    i=0, j=0, k=0, ix, iy, iz;
  /* 
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

			if(rk==0.0){
              rky[i] = 0.0;
              rkx[i] = 0.0;
              rkz[i] = 0.0;
			}else{
              rky[i] = ky/rk;
              rkx[i] = kx/rk;
              rkz[i] = kz/rk;
			}
            i++;
         }
      }
   }
   */

   for(iy=0; iy < nky; iy++)
   {
     ky = ky0+iy*dky;
     for(ix=0; ix < nkx; ix++)
     {
       kx = kx0+ix*dkx;
       for (iz=0; iz < nkz; iz++)
       {
            kz = kz0+iz*dkz;
			if(ky==0.0&&kx==0.0&&kz==0.0){
              rky[i] = 0.0;
              rkx[i] = 0.0;
              rkz[i] = 0.0;
			}else{
              k2 = ky*ky+kx*kx+kz*kz;
              rk = sqrt(k2);
              rky[i] = ky/rk;
              rkx[i] = kx/rk;
              rkz[i] = kz/rk;
			}
            i++;
         }
      }
   }

   t2=clock();
   timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
   sf_warning("CPU time for prereparing for div-curl-operation: %f(second)",timespent);

   /****************begin to input wavefield (snapshot) and separation ****************/
   /****************begin to input wavefield (snapshot) and separation ****************/
   int *ijky = sf_intalloc(nky);
   int *ijkx = sf_intalloc(nkx);
   int *ijkz = sf_intalloc(nkz);

   ikxikyikz(ijkx, ijky, ijkz, nkx, nky, nkz);

    float *px, *py, *pz,*p,*pp;

    px=sf_floatalloc(nxyz);
    py=sf_floatalloc(nxyz);
    pz=sf_floatalloc(nxyz);
    p=sf_floatalloc(nxyz);
    pp=sf_floatalloc(nxyz);

    int iflag=0;

	sf_floatread(px, nxyz, Fx);
	sf_floatread(py, nxyz, Fy);
	sf_floatread(pz, nxyz, Fz);

	sf_file Fp, Fsv, Fsh;
    Fp = sf_output("out");
    Fsv = sf_output("ElasticSepSV");
    Fsh = sf_output("ElasticSepSH");

	puthead3x(Fp,  nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
	puthead3x(Fsv, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);
	puthead3x(Fsh, nz, nx, ny, dz/1000, dx/1000, dy/1000, fz/1000, fx/1000, fy/1000);

    sf_warning("separate qP-wave based on div-curl-operation."); 
    // separate qP wave  
    for(k=0;k<nxyz;k++){
		p[k] = px[k];
	}
    sepdiv3dD(rkx,p,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] = p[k];

    for(k=0;k<nxyz;k++)
		p[k] = py[k];
    sepdiv3dD(rky,p,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] += p[k];

    for(k=0;k<nxyz;k++)
		p[k] = pz[k];
    sepdiv3dD(rkz,p,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] += p[k];

	sf_floatwrite(pp, nxyz, Fp);
    FILE *fp;
    fp=fopen("snap3dep.pc","wb");
    fwrite(pp,sizeof(float),nxyz,fp);
    fclose(fp);

   double *rksx, *rksy, *rksz, ss;
   rksx=(double*)calloc(sizeof(double),nxyz);
   rksy=(double*)calloc(sizeof(double),nxyz);
   rksz=(double*)calloc(sizeof(double),nxyz);

    // separate SH wave  
    sf_warning("separate SH-wave based on div-curl-operation."); 
    for(k=0;k<nxyz;k++){
		double r = rky[k]*rky[k]+rkx[k]*rkx[k];
		if(r==0)
		   ss = 0.0;
		else
		   ss = sqrt(1.0-rkz[k]*rkz[k])/sqrt(r);
		rksy[k]=-rky[k]*ss;
		rksx[k]= rkx[k]*ss;

		p[k] = px[k];
	}

    sepdiv3dD(rksy,p,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] = p[k];

    for(k=0;k<nxyz;k++) p[k] = py[k];
    sepdiv3dD(rksx,p,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] += p[k];

	sf_floatwrite(pp, nxyz, Fsh);
    fp=fopen("snap3desh.pc","wb");
    fwrite(pp,sizeof(float),nxyz,fp);
    fclose(fp);

    // separate SV wave  
    sf_warning("separate SV-wave based on div-curl-operation."); 
    for(k=0;k<nxyz;k++){
		rksx[k]= -rkx[k]*rkz[k];
		rksy[k]= -rky[k]*rkz[k];
		rksz[k]= rky[k]*rky[k]+rkx[k]*rkx[k];
		double r = rksx[k]*rksx[k]+rksy[k]*rksy[k]+rksz[k]*rksz[k];
		if(r==0.0)
		   ss = 0.0;
		else
		   ss = sqrt(1.0-rkz[k]*rkz[k])/sqrt(r);
		rksx[k] *= ss;
		rksy[k] *= ss;
		rksz[k] *= ss;
	}
    sepdiv3dD(rksx,px,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] = px[k];

    sepdiv3dD(rksy,py,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] += py[k];

    sepdiv3dD(rksz,pz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,iflag);
    for(k=0;k<nxyz;k++) pp[k] += pz[k];

	sf_floatwrite(pp, nxyz, Fsv);
    fp=fopen("snap3desv.pc","wb");
    fwrite(pp,sizeof(float),nxyz,fp);
    fclose(fp);

    free(pp);
    free(p);
    free(px);
    free(py);
    free(pz);
    free(rkx);
    free(rky);
    free(rkz);
    free(rksx);
    free(rksy);
    free(rksz);

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wave-modes separation.: %f(second)",timespent);

    sf_warning("-------sucessful ending --------");
    exit(0);
}
