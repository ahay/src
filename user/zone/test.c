/* Forward propagating using original elastic equation of displacement in ORT media */ 
/* *************************************************************************
 * Copyright: Tongji University (Jiubing Cheng and Tengfei Wang)
 *     2012.3.2
 * *************************************************************************/
/*
  Copyright (C) 2012 Tongji University, Shanghai, China 
  Authors: Jiubing Cheng and Tengfei Wang
  Modified: Yanadet Sripanich (The University of Texas at Austin) 
     
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "_cjb.h"
#include "_fd.h"

#include "zero.h"


void fwportelastic3d(float dt2,float***p1,float***p2,float***p3,float***q1,float***q2,float***q3,
		     float***r1,float***r2,float***r3,float***px_tmp,float ***pz_tmp,
		     float***qx_tmp,float***qz_tmp,
		     float***rx_tmp,float***rz_tmp,
                     float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,
		     float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                     float dx, float dy, float dz, int nxpad, int nypad, int nzpad, 
                     float ***c11, float ***c22, float ***c33, float ***c12, float ***c13, float ***c23, float ***c44, float ***c55, float ***c66,
		     float ***phaix, float ***phaiy, float ***phaiz)
/*< fwpttielastic3d: forward-propagating using original elastic equation of displacement in 3D TTI media>*/
{
    int   i,j,k,l;
    float pxy=0, pxz=0, pyz=0, qxy, qxz, qyz, rxy, rxz, ryz;
    float phix,phiy,phiz;
    float sinphiz,cosphiz,sinphiy,cosphiy,sinphix,cosphix;
    float r11, r12, r13, r21, r22, r23, r31, r32, r33;
    float a11, a22, a33, a12, a13, a23, a44, a55, a66;
    float hpy, hqy, hry, hpx, hqx, hrx, hpz, hqz, hrz;

    float px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2;
    float qxy1, rxz1, pxy1, ryz1, pxz1, qyz1;
		
    zero3float(px_tmp,nzpad,nxpad,nypad);	
    zero3float(pz_tmp,nzpad,nxpad,nypad);	
    zero3float(qx_tmp,nzpad,nxpad,nypad);	
    zero3float(qz_tmp,nzpad,nxpad,nypad);	
    zero3float(rx_tmp,nzpad,nxpad,nypad);	
    zero3float(rz_tmp,nzpad,nxpad,nypad);	

#ifdef _OPENMP
#pragma omp parallel for private(i,j,k,l)				\
    schedule(dynamic)							\
    shared(p2,q2,r2,							\
	   px_tmp,pz_tmp,						\
	   qx_tmp,qz_tmp,						\
	   rx_tmp,rz_tmp,						\
	   coeff_1dx,coeff_1dy,coeff_1dz,nxpad,nypad,nzpad,dx,dy,dz)
#endif
    for(k=0;k<nypad;k++)
	for(i=0;i<nxpad;i++)
	    for(j=0;j<nzpad;j++)
		for(l=-_mix;l<=_mix;l++)
		{
		    if(i+l>=0&&i+l<nxpad){
			px_tmp[k][i][j]+=coeff_1dx[l+_mix]*p2[k][i+l][j]/2.0/dx;
			qx_tmp[k][i][j]+=coeff_1dx[l+_mix]*q2[k][i+l][j]/2.0/dx;
			rx_tmp[k][i][j]+=coeff_1dx[l+_mix]*r2[k][i+l][j]/2.0/dx;
		    }
		    if(j+l>=0&&j+l<nzpad){
			pz_tmp[k][i][j]+=coeff_1dz[l+_mix]*p2[k][i][j+l]/2.0/dz;
			qz_tmp[k][i][j]+=coeff_1dz[l+_mix]*q2[k][i][j+l]/2.0/dz;
			rz_tmp[k][i][j]+=coeff_1dz[l+_mix]*r2[k][i][j+l]/2.0/dz;
		    }
		}

#ifdef _OPENMP
#pragma omp parallel for private(i,j,k,l,hpy, hqy, hry, hpx, hqx, hrx, hpz, hqz, hrz, px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2, qxy1, rxz1, pxy1, ryz1, pxz1, qyz1) \
    schedule(dynamic)							\
    shared(p1,p2,p3,q1,q2,q3,r1,r2,r3,					\
	   px_tmp,pz_tmp,						\
	   qx_tmp,qz_tmp,						\
	   rx_tmp,rz_tmp,						\
	   coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz,	\
	   c11,c22,c33,c12,c13,c23,c44,c55,c66,phaix,phaiy,phaiz,dt2)
#endif
    for(k=0;k<nypad;k++)
	for(i=0;i<nxpad;i++)
	    for(j=0;j<nzpad;j++)
	    {
		
		a11=c11[k][i][j];
		a22=c22[k][i][j];
		a33=c33[k][i][j];
		a12=c12[k][i][j];
		a13=c13[k][i][j];
		a23=c23[k][i][j];
		a44=c44[k][i][j];
		a55=c55[k][i][j];
		a66=c66[k][i][j];
		phix=phaix[k][i][j];
		phiy=phaiy[k][i][j];
		phiz=phaiz[k][i][j];

		cosphix=cos(phix);
		sinphix=sin(phix);
		cosphiy=cos(phiy);
		sinphiy=sin(phiy);
		cosphiz=cos(phiz);
		sinphiz=sin(phiz);

		r11=cosphiy*cosphiz;
		r12=cosphix*sinphiz + sinphix*sinphiy*cosphiz;
		r13=sinphix*sinphiz - cosphix*sinphiy*cosphiz;
		r21=-cosphiy*sinphiz;
		r22=cosphix*cosphiz - sinphix*sinphiy*sinphiz;
		r23=sinphix*cosphiz + cosphix*sinphiy*sinphiz;
		r31=sinphiy;
		r32=-sinphix*cosphiy;
		r33=cosphix*cosphiy;


		pxy=0;
		qxy=0;
		rxy=0;
		pxz=0;
		qxz=0;
		rxz=0;
		pyz=0;
		qyz=0;
		ryz=0;
		for(l=-_mix;l<=_mix;l++)
		{
		    if(k+l>=0&&k+l<nypad)
		    {
			pxy+=coeff_1dy[l+_mix]*px_tmp[k+l][i][j]/2.0/dy;
			pyz+=coeff_1dy[l+_mix]*pz_tmp[k+l][i][j]/2.0/dy;

			qxy+=coeff_1dy[l+_mix]*qx_tmp[k+l][i][j]/2.0/dy;
			qyz+=coeff_1dy[l+_mix]*qz_tmp[k+l][i][j]/2.0/dy;

			rxy+=coeff_1dy[l+_mix]*rx_tmp[k+l][i][j]/2.0/dy;
			ryz+=coeff_1dy[l+_mix]*rz_tmp[k+l][i][j]/2.0/dy;
		    }
		    if(j+l>=0&&j+l<nzpad)
		    {
			pxz+=coeff_1dz[l+_mix]*px_tmp[k][i][j+l]/2.0/dz;
			qxz+=coeff_1dz[l+_mix]*qx_tmp[k][i][j+l]/2.0/dz;
			rxz+=coeff_1dz[l+_mix]*rx_tmp[k][i][j+l]/2.0/dz;
		    }
		}
		/* 2nd-order derivatives */
		hpy =0;
		hqy =0;
		hry =0;
		hpx =0;
		hqx =0;
		hrx =0;
		hpz =0;
		hqz =0;
		hrz =0;
		for(l=-_m;l<=_m;l++)
		{
		    if(k+l>=0&&k+l<nypad)
		    {
			hpy +=coeff_2dy[l+_m]*p2[k+l][i][j];
			hqy +=coeff_2dy[l+_m]*q2[k+l][i][j];
			hry +=coeff_2dy[l+_m]*r2[k+l][i][j];
		    }
		    if(i+l>=0&&i+l<nxpad)
		    {
			hpx +=coeff_2dx[l+_m]*p2[k][i+l][j];
			hqx +=coeff_2dx[l+_m]*q2[k][i+l][j];
			hrx +=coeff_2dx[l+_m]*r2[k][i+l][j];
		    }
		    if(j+l>=0&&j+l<nzpad)
		    {
			hpz +=coeff_2dz[l+_m]*p2[k][i][j+l];
			hqz +=coeff_2dz[l+_m]*q2[k][i][j+l];
			hrz +=coeff_2dz[l+_m]*r2[k][i][j+l];
		    }
		}

		/* x-component */
		px2 = r11*r11*hpx+r21*r21*hpy+r31*r31*hpz+2*r11*r21*pxy+2*r11*r31*pxz+2*r21*r31*pyz;
		py2 = r12*r12*hpx+r22*r22*hpy+r32*r32*hpz+2*r12*r22*pxy+2*r12*r32*pxz+2*r22*r32*pyz;
		pz2 = r13*r13*hpx+r23*r23*hpy+r33*r33*hpz+2*r13*r23*pxy+2*r13*r33*pxz+2*r23*r33*pyz;

		qxy1 = r11*r12*hqx+r21*r22*hqy+r31*r32*hqz+(r11*r22+r21*r12)*qxy+(r11*r32+r31*r12)*qxz+(r21*r32+r31*r22)*qyz;
		rxz1 = r11*r13*hrx+r21*r23*hry+r31*r33*hrz+(r11*r23+r21*r13)*rxy+(r11*r33+r31*r13)*rxz+(r21*r33+r31*r23)*ryz;

		p3[k][i][j]=2*p2[k][i][j] - p1[k][i][j] + dt2*( a11*px2 +  a66*py2 + a55*pz2 + (a12+a66)*qxy1 + (a13+a55)*rxz1);

		/* y-component */
		qx2 = r11*r11*hqx+r21*r21*hqy+r31*r31*hqz+2*r11*r21*qxy+2*r11*r31*qxz+2*r21*r31*qyz;
		qy2 = r12*r12*hqx+r22*r22*hqy+r32*r32*hqz+2*r12*r22*qxy+2*r12*r32*qxz+2*r22*r32*qyz;
		qz2 = r13*r13*hqx+r23*r23*hqy+r33*r33*hqz+2*r13*r23*qxy+2*r13*r33*qxz+2*r23*r33*qyz;

		pxy1 = r11*r12*hpx+r21*r22*hpy+r31*r32*hpz+(r11*r22+r21*r12)*pxy+(r11*r32+r31*r12)*pxz+(r21*r32+r31*r22)*pyz;
		ryz1 = r12*r13*hrx+r22*r23*hry+r32*r33*hrz+(r12*r23+r22*r13)*rxy+(r12*r33+r32*r13)*rxz+(r22*r33+r32*r23)*ryz;

		q3[k][i][j]=2*q2[k][i][j] - q1[k][i][j] + dt2*( a66*qx2 +  a22*qy2 + a44*qz2 + (a12+a66)*pxy1 + (a23+a44)*ryz1);

		/* z-component */
		rx2 = r11*r11*hrx+r21*r21*hry+r31*r31*hrz+2*r11*r21*rxy+2*r11*r31*rxz+2*r21*r31*ryz;
		ry2 = r12*r12*hrx+r22*r22*hry+r32*r32*hrz+2*r12*r22*rxy+2*r12*r32*rxz+2*r22*r32*ryz;
		rz2 = r13*r13*hrx+r23*r23*hry+r33*r33*hrz+2*r13*r23*rxy+2*r13*r33*rxz+2*r23*r33*ryz;

		pxz1 = r11*r13*hpx+r21*r23*hpy+r31*r33*hpz+(r11*r23+r21*r13)*pxy+(r11*r33+r31*r13)*pxz+(r21*r33+r31*r23)*pyz;
		qyz1 = r12*r13*hqx+r22*r23*hqy+r32*r33*hqz+(r12*r23+r22*r13)*qxy+(r12*r33+r32*r13)*qxz+(r22*r33+r32*r23)*qyz;

		r3[k][i][j]=2*r2[k][i][j] - r1[k][i][j] + dt2*( a55*rx2 +  a44*ry2 + a33*rz2 + (a13+a55)*pxz1 + (a23+a44)*qyz1);
	    }
}

