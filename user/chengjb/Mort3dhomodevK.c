/* 3D three-components projection deviation correction operators calculation in
 * homogeneous orthorhombic media
   Copyright (C) 2012 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng and Tengfei Wang 
     
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

#include <rsf.h>
#include <assert.h>

/* Joe's taper */
#define TAPER(k) (0.5*(1+cosf(k))) 

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from C programs */
#include "zero.h"
#include "kykxkztaper.h"
#include "puthead.h"
#include "kykxkz2yxz.h"
#include "eigen3x3.h"

int main(int  argc,char **argv)
{
    sf_init(argc,argv);

    float   vp0,vs0,ga1,ep1,de1,ga2,ep2,de2,de3,alpha,the,phi;
    float   dx,dy,dz,dkx,dky,dkz,fkx,fky,fkz;

    int     i,j,k,ix,iy,hny,hnx,hnz,ny,nx,nz,itaper;

    double  sinx, siny, sinz, sinx2, siny2, sinz2, k2, vp2, vs2, kx, ky, kz;

    clock_t t1, t2;
    float   timespent;

    t1=clock();

    /* time samping paramter */
    if (!sf_getfloat("vp0",&vp0)) vp0=3000.0;
    if (!sf_getfloat("vs0",&vs0)) vs0=1500.0;
    if (!sf_getfloat("de1",&de1)) de1=0.05;
    if (!sf_getfloat("de2",&de2)) de2=-0.05;
    if (!sf_getfloat("de3",&de3)) de3=0.15;
    if (!sf_getfloat("ep1",&ep1)) ep1=0.2;
    if (!sf_getfloat("ep2",&ep2)) ep2=0.05;
    if (!sf_getfloat("ga1",&ga1)) ga1=0.1;
    if (!sf_getfloat("ga2",&ga2)) ga2=0.1;

    if (!sf_getfloat("alpha",&alpha)) alpha=0.;
    if (!sf_getfloat("the",&the)) the=0.;
    if (!sf_getfloat("phi",&phi)) phi=0.;

    if (!sf_getint("hnx",&hnx)) hnx=250;
    if (!sf_getint("hny",&hny)) hny=250;
    if (!sf_getint("hnz",&hnz)) hnz=250;

    if (!sf_getfloat("dx",&dx)) dx=10.;
    if (!sf_getfloat("dy",&dy)) dy=10.;
    if (!sf_getfloat("dz",&dz)) dz=10.;

    if (!sf_getint("itaper",&itaper)) itaper=1;
    nx=2*hnx+1;
    ny=2*hny+1;
    nz=2*hnz+1;

    sf_warning("itaper=%d ",itaper);
    sf_warning("nx=%d ",nx);
    sf_warning("ny=%d ",ny);
    sf_warning("nz=%d ",nz);
    sf_warning("dx=%f ",dx);
    sf_warning("dy=%f ",dy);
    sf_warning("dz=%f ",dz);

    dkx=2*SF_PI/dx/nx;
    dkz=2*SF_PI/dz/nz;
    dky=2*SF_PI/dy/ny;

    fkx=-SF_PI/dx;
    fky=-SF_PI/dy;
    fkz=-SF_PI/dz;

    sf_warning("vp0=%f ",vp0);
    sf_warning("vs0=%f ",vs0);
    sf_warning("de1=%f ",de1);
    sf_warning("de2=%f ",de2);
    sf_warning("de3=%f ",de3);
    sf_warning("ep1=%f ",ep1);
    sf_warning("ep2=%f ",ep2);
    sf_warning("ga1=%f ",ga1);
    sf_warning("ga2=%f ",ga2);

    sf_warning("alpha=%f ",alpha);
    sf_warning("the=%f ",the);
    sf_warning("phi=%f ",phi);

    alpha *= SF_PI/180.0;
    the *= SF_PI/180.0;
    phi *= SF_PI/180.0;
    
    char    jobz='V';  /* for SVD */
    char    uplo='U';  /* for SVD */
    int     M=3;  /* for SVD */
    int     LDA=M;  /* for SVD */
    int     LWORK=8*M;  /* for SVD */
    int     INFO;  /* for SVD */
    double  *Chr, *ww, *work;  /* Lapack SVD array */

    //Chr=sf_doublealloc(M*M);    // no code in Madagascar
    //ww=sf_doublealloc(M*M);    // no code in Madagascar
    //work=sf_doublealloc(LWORK);    // no code in Madagascar

    Chr=calloc(sizeof(double),M*M);
    ww=calloc(sizeof(double),M*M);
    work=calloc(sizeof(double),LWORK);
    
    //double A[3][3],w[3],Q[3][3];

    float*** apvx=sf_floatalloc3(nz,nx,ny);
    float*** apvy=sf_floatalloc3(nz,nx,ny);
    float*** apvz=sf_floatalloc3(nz,nx,ny);

    vp2=vp0*vp0;
    vs2=vs0*vs0;
    ep1=1+2*ep1;
    de1=1+2*de1;
    ga1=1+2*ga1;
    ep2=1+2*ep2;
    de2=1+2*de2;
    ga2=1+2*ga2;
    de3=1+2*de3;

    double a11, a22, a33, a44, a55, a66, a12a66, a23a44, a13a55;

    a11 = ep2*vp2;
    a22 = ep1*vp2;
    a33 = vp2;
    a44 = ga1/ga2*vs2;
    a55 = vs2;
    a66 = ga1*vs2;
    a12a66 = sqrt((a11-a66)*(de3*a11-a66));//a12+a66
    a13a55 = sqrt((a33-a55)*(de2*a33-a55));//a13+a55
    a23a44 = sqrt((a33-a44)*(de1*a33-a44));//a23+a44

    int jhny, ihnx, khnz;

    zero3float(apvy, nz, nx, ny);
    zero3float(apvx, nz, nx, ny);
    zero3float(apvz, nz, nx, ny);

    int      iz;
    float    rkx, rky, rkz;

    float*** taper=sf_floatalloc3(nz,nx,ny);
    
    for( j=-hny; j<=hny; j++){
          jhny=j+hny;
	  ky=j*dky;
          rky=2*SF_PI*j/ny;
         //sf_warning("j=%d jhny=%d ky=%f ",j,jhny,ky);
	  for( i=-hnx; i<=hnx; i++ ){
                ihnx=i+hnx;
		kx=i*dkx;
                rkx=2*SF_PI*i/nx;
		for( k=-hnz; k<=hnz; k++)
		{
                        khnz=k+hnz;
			kz=k*dkz;
                        rkz=2*SF_PI*k/nz;

                        taper[jhny][ihnx][khnz]=pow((TAPER(rky)*TAPER(rkx)*TAPER(rkz)), 1.0/100.0);

			if(i==0 && j==0 && k==0)
                        {
			   apvy[jhny][ihnx][khnz]=1.0;
			   apvx[jhny][ihnx][khnz]=1.0;
			   apvz[jhny][ihnx][khnz]=1.0;
			   continue;
                        }

                        //sf_warning("j= %d i= %d  k=%d ",j,i,k);

                         if(kx==0.0)   //move this sample from zero for cintinuous spectrum
                             kx = 0.0000000001*dkx;
                         if(ky==0.0)   //move this sample from zero for cintinuous spectrum
                             ky = 0.0000000001*dky;
                         if(kz==0.0)   //move this sample from zero for cintinuous spectrum
                             kz = 0.0000000001*dkz;

			k2=sqrt(kx*kx+ky*ky+kz*kz);
			sinx=kx/k2;
			siny=ky/k2;
			sinz=kz/k2;

			sinx2=sinx*sinx;
			siny2=siny*siny;
			sinz2=sinz*sinz;

                        double upy, upx, upz;
/*
			A[0][0] = a11*sinx2 + a66*siny2 + a55*sinz2;
			A[0][1] = A[0][1] = a12a66*sinx*siny;
			A[0][2] = A[2][0] = a13a55*sinx*sinz;
			A[1][1] = a66*sinx2 + a22*siny2 + a44*sinz2;
			A[1][2] = A[2][1] = a23a44*siny*sinz;
			A[2][2] = a55*sinx2 + a44*siny2 + a33*sinz2;

                        // Hybrid algorithm of eigeinvalue solution for hermitian 3X3 matrices
                        // Faster than LAPACK's ssyev routine about 7 times
                        // But less accurate !!!
                        if( dsyevh3(A, Q, w) != 0) {
                          sf_warning("dsyevh3 error !!!");
                          exit(0);
                        }

                        upx=Q[0][0];
                        upy=Q[1][0];
                        upz=Q[2][0];

                        // get the closest direction to k 
                        if(upx*sinx + upy*siny+ upz*sinz < 0.) {
                            upx=-Q[0][0];
                            upy=-Q[1][0];
                            upz=-Q[2][0];
                        }
                        sf_warning("v1=%f v2=%f v3=%f ",sqrt(w[0]), sqrt(w[1]),sqrt(w[2]));
                        sf_warning("eigeinvector1: %f %f %f",(float)Q[0][0],(float)Q[1][0],(float)Q[2][0]);
                        sf_warning("eigeinvector2: %f %f %f",(float)Q[0][1],(float)Q[1][1],(float)Q[2][1]);
                        sf_warning("eigeinvector3: %f %f %f",(float)Q[0][2],(float)Q[1][2],(float)Q[2][2]);
*/
			Chr[0] = a11*sinx2 + a66*siny2 + a55*sinz2;
			Chr[1] = a12a66*sinx*siny;
			Chr[2] = a13a55*sinx*sinz;
			Chr[3] = Chr[1];
			Chr[4] = a66*sinx2 + a22*siny2 + a44*sinz2;
			Chr[5] = a23a44*siny*sinz;
			Chr[6] = Chr[2];
			Chr[7] = Chr[5];
			Chr[8] = a55*sinx2 + a44*siny2 + a33*sinz2;

                        // LAPACK's ssyev routine (slow but accurate)
                        // double type is more accurate than single type
			dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);

                        upx=Chr[6];
                        upy=Chr[7];
                        upz=Chr[8];
                        if(upx*sinx + upy*siny+ upz*sinz < 0.) {
                            upx=-Chr[6];
                            upy=-Chr[7];
                            upz=-Chr[8];
                        }
/*
                        sf_warning("v1=%f v2=%f v3=%f ",sqrt(ww[0]), sqrt(ww[1]),sqrt(ww[2]));
                        sf_warning("eigeinvector1: %f %f %f",Chr[0],Chr[1],Chr[2]);
                        sf_warning("eigeinvector2: %f %f %f",Chr[3],Chr[4],Chr[5]);
                        sf_warning("eigeinvector3: %f %f %f",Chr[6],Chr[7],Chr[8]);
                        sf_warning("=======================================================");
			apvy[jhny][ihnx][khnz]=(float)fabs(upy/siny);
			apvx[jhny][ihnx][khnz]=(float)fabs(upx/sinx);
			apvz[jhny][ihnx][khnz]=(float)fabs(upz/sinz);
*/
			apvy[jhny][ihnx][khnz]=(float)(upy/siny);
			apvx[jhny][ihnx][khnz]=(float)(upx/sinx);
			apvz[jhny][ihnx][khnz]=(float)(upz/sinz);

		}// k loop
           } // i loop
      } // j loop

/* Wang Tengfei
 */
    for( j=0; j<ny ; j++)
         for( i=0; i<nx; i++ )
            for( k=0; k<nz ; k++)
            {
                if(fabs(apvy[j][i][k])==0)
                {
                    if(j==0)
                       apvy[j][i][k]=apvy[j+1][i][k];
                    else if(j==ny-1)
                       apvy[j][i][k]=apvy[j-1][i][k];
                    else
                       apvy[j][i][k]=(apvy[j+1][i][k] + apvy[j-1][i][k])/2.0;
                }
                if(fabs(apvx[j][i][k])==0)
                {
                    if(i==0)
                       apvx[j][i][k]=apvx[j][i+1][k];
                    else if(i==nx-1)
                       apvx[j][i][k]=apvx[j][i-1][k];
                    else
                       apvx[j][i][k]=(apvx[j][i+1][k] + apvx[j][i-1][k])/2.0;
                }
                if(fabs(apvz[j][i][k])==0)
                {
                    if(k==0)
                       apvz[j][i][k]=apvz[j][i][k+1] ;
                    else if(k==nz-1)
                       apvz[j][i][k]=apvz[j][i][k-1] ;
                    else
                       apvz[j][i][k]=(apvz[j][i][k+1] + apvz[j][i][k-1])/2.0;
                }
            }

//  Cheng Jiubing
    for( j=0; j<ny ; j++)
         for( i=0; i<nx; i++ )
         {
             apvy[j][i][hnz]=(apvy[j][i][hnz+1] + apvy[j][i][hnz-1])/2.0;
             apvx[j][i][hnz]=(apvx[j][i][hnz+1] + apvx[j][i][hnz-1])/2.0;
             apvz[j][i][hnz]=(apvz[j][i][hnz+1] + apvz[j][i][hnz-1])/2.0;
         }
    for( j=0; j<ny ; j++)
         for( k=0; k<nz; k++ )
         {
             apvy[j][hnx][k]=(apvy[j][hnx+1][k] + apvy[j][hnx-1][k])/2.0;
             apvx[j][hnx][k]=(apvx[j][hnx+1][k] + apvx[j][hnx-1][k])/2.0;
             apvz[j][hnx][k]=(apvz[j][hnx+1][k] + apvz[j][hnx-1][k])/2.0;
         }
    for( i=0; i<nx ; i++)
         for( k=0; k<nz; k++ )
         {
             apvy[hny][i][k]=(apvy[hny+1][i][k] + apvy[hny-1][i][k])/2.0;
             apvx[hny][i][k]=(apvx[hny+1][i][k] + apvx[hny-1][i][k])/2.0;
             apvz[hny][i][k]=(apvz[hny+1][i][k] + apvz[hny-1][i][k])/2.0;
         }
//
    sf_file Fo1, Fo2, Fo3, Fo4;

    Fo1 = sf_output("out");    /* deviation operator's x-component in (ky,kx,kz) domain */
    Fo2 = sf_output("apvy");   /* deviation operator's y-component in (ky,kx,kz) domain */
    Fo3 = sf_output("apvz");   /* deviation operator's z-component in (ky,kx,kz) domain */
    Fo4 = sf_output("taper");  /* deviation operator's z-component in (ky,kx,kz) domain */

    puthead3kx(Fo1, nz, nx, ny, dkz, dkx, dky, fkz, fkx, fky);
    puthead3kx(Fo2, nz, nx, ny, dkz, dkx, dky, fkz, fkx, fky);
    puthead3kx(Fo3, nz, nx, ny, dkz, dkx, dky, fkz, fkx, fky);
    puthead3kx(Fo4, nz, nx, ny, dkz, dkx, dky, fkz, fkx, fky);

    if(itaper==1)
    {
      for(iy=0;iy<ny;iy++)
      for(ix=0;ix<nx;ix++)
      for(iz=0;iz<nz;iz++)
      {
       apvx[iy][ix][iz] *= taper[iy][ix][iz];
       apvy[iy][ix][iz] *= taper[iy][ix][iz];
       apvz[iy][ix][iz] *= taper[iy][ix][iz];
      }
    }

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatwrite(apvx[iy][ix],nz,Fo1);
      sf_floatwrite(apvy[iy][ix],nz,Fo2);
      sf_floatwrite(apvz[iy][ix],nz,Fo3);
      sf_floatwrite(taper[iy][ix],nz,Fo4);
    }
    free(**taper);
    free(**apvx);
    free(**apvy);
    free(**apvz);

    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("CPU time for 3D ORT deviation operators: %f(second)",timespent);

   /****************End of Calculating Polarization Vector****************/
    free(Chr);
    free(ww);
    free(work);

return 0;
}

