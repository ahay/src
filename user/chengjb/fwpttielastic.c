/*************************************************************************
 * * Forward propagating using original elastic equation of displacement 
 *   in VTI media
 * 
 *     Copyright: Tongji University (Jiubing Cheng and Tengfei Wang)
 *     2012.3.2
 * *************************************************************************/
/*
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
#include <omp.h>

#include "_cjb.h"
#include "_fd.h"

#include "zero.h"

void fwpttielastic(float dt2, float** p1,float** p2,float** p3, float** q1,float** q2,float** q3,
                   float* coeff_2dx,float* coeff_2dz, float* coeff_1dx,float* coeff_1dz,
                   float dx, float dz, int nx, int nz, int nxpad, int nzpad, 
                   float **vp0,float **vs0, float **epsilon,float **delta, float **theta)
/*< fwpttielastic: forward-propagating using original elastic equation of displacement in 2D TTI media>*/
{
        int   i,j,l;
		
	float **px_tmp=sf_floatalloc2(nzpad,nxpad);
	float **qx_tmp=sf_floatalloc2(nzpad,nxpad);

	zero2float(px_tmp,nzpad,nxpad);	
	zero2float(qx_tmp,nzpad,nxpad);	

#pragma omp parallel for private(i,j,l) \
	    schedule(dynamic) \
        shared(p2,q2,px_tmp,qx_tmp,coeff_1dx,dx)
    for(i=_m;i<nx+_m;i++)
	for(j=_m;j<nz+_m;j++)
	{
		for(l=-_mix;l<=_mix;l++)
		{
            int lm=l+_mix;
			px_tmp[i][j]+=coeff_1dx[lm]*p2[i+l][j]/2.0/dx;
			qx_tmp[i][j]+=coeff_1dx[lm]*q2[i+l][j]/2.0/dx;
		}
	}

#pragma omp parallel for private(i,j,l) \
	    schedule(dynamic) \
	    shared(p1,p2,p3,q1,q2,q3,\
		px_tmp,qx_tmp,\
		coeff_1dx,coeff_1dz,coeff_2dx,coeff_2dz,\
	    vp0,vs0,epsilon,delta,theta,dt2)
        for(i=_m;i<nx+_m;i++)
        {
           int im=i-_m;
	   for(j=_m;j<nz+_m;j++)
	   {
               int jm=j-_m;

               float vp2,vs2,ep,de,vpx2,vpn2,coef, the;
               float sinthe,costhe,cos2,sin2,sin2a,cos_sin;
               float px,pxz,qxz,qx,px1, qxz1, qx1, pxz1,hpx,hqx,hpz,hqz;

               vp2=vp0[im][jm]*vp0[im][jm];
               vs2=vs0[im][jm]*vs0[im][jm];
               ep=1+2*epsilon[im][jm];
               de=1+2*delta[im][jm];
               the=theta[im][jm];

               costhe=cos(the);
               sinthe=sin(the);
               cos2=costhe*costhe;
               sin2=sinthe*sinthe;
               cos_sin=costhe*sinthe;
               sin2a=2*cos_sin;

	       vpx2=vp2*ep;
	       vpn2=vp2*de;
               coef=sqrt((vp2-vs2)*(vpn2-vs2));

		pxz=0;
		qxz=0;
                for(l=-_mix;l<=_mix;l++)
                {
                     int lm=l+_mix;
                     pxz+=coeff_1dz[lm]*px_tmp[i][j+l]/2.0/dz;
                     qxz+=coeff_1dz[lm]*qx_tmp[i][j+l]/2.0/dz;
                }

		hpx=0;
                hpz=0;
		hqx=0;
                hqz=0;
		for(l=-_m;l<=_m;l++)
		{
                     int lm=l+_m;
                     hpx+=coeff_2dx[lm]*p2[i+l][j];
                     hqx+=coeff_2dx[lm]*q2[i+l][j];
                     hpz+=coeff_2dz[lm]*p2[i][j+l];
                     hqz+=coeff_2dz[lm]*q2[i][j+l];
		}

                px  = cos2*hpx + sin2*hpz + sin2a*pxz;
                px1 = sin2*hpx + cos2*hpz - sin2a*pxz;
                qxz1 = -cos_sin*hqx + cos_sin*hqz + (cos2-sin2)*qxz;

                p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*( vpx2*px +  vs2*px1 + sqrt((vp2-vs2)*(vpn2-vs2))*qxz1);

                qx  = cos2*hqx + sin2*hqz + sin2a*qxz;
                qx1 = sin2*hqx + cos2*hqz - sin2a*qxz;
                pxz1 = -cos_sin*hpx + cos_sin*hpz + (cos2-sin2)*pxz;

                q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*( vs2*qx +  vp2*qx1 + sqrt((vp2-vs2)*(vpn2-vs2))*pxz1);
          }
	}

	free(*px_tmp);	
	free(*qx_tmp);	
}

void fwpttielastic3d(float dt2,float***p1,float***p2,float***p3,float***q1,float***q2,float***q3,
		             float***r1,float***r2,float***r3,float***px_tmp,float ***pz_tmp,
					 float***qx_tmp,float***qz_tmp,
					 float***rx_tmp,float***rz_tmp,
                     float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,
					 float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                     float dx, float dy, float dz, int nxpad, int nypad, int nzpad, 
                     float ***vp0,float ***vs0, float ***epsilon, float ***delta, float ***gama,
					 float ***theta, float ***phai)
/*< fwpttielastic3d: forward-propagating using original elastic equation of displacement in 3D TTI media>*/
{
    int   i,j,k,l;
		
	zero3float(px_tmp,nzpad,nxpad,nypad);	
	zero3float(pz_tmp,nzpad,nxpad,nypad);	
	zero3float(qx_tmp,nzpad,nxpad,nypad);	
	zero3float(qz_tmp,nzpad,nxpad,nypad);	
	zero3float(rx_tmp,nzpad,nxpad,nypad);	
	zero3float(rz_tmp,nzpad,nxpad,nypad);	

#pragma omp parallel for private(i,j,k,l) \
	    schedule(dynamic) \
        shared(p2,q2,r2, \
		px_tmp,pz_tmp, \
	    qx_tmp,qz_tmp,\
		rx_tmp,rz_tmp, \
		coeff_1dx,coeff_1dy,coeff_1dz,nxpad,nypad,nzpad,dx,dy,dz)
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

#pragma omp parallel for private(i,j,k,l) \
	    schedule(dynamic) \
	    shared(p1,p2,p3,q1,q2,q3,r1,r2,r3,\
		px_tmp,pz_tmp, \
	    qx_tmp,qz_tmp,\
		rx_tmp,rz_tmp, \
		coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz,\
	    vp0,vs0,epsilon,delta,gama,theta,phai,dt2)
	for(k=0;k<nypad;k++)
	for(i=0;i<nxpad;i++)
		for(j=0;j<nzpad;j++)
		{
			float vp2,vs2,ep,de,ga,vpn2,the,phi;
			float sinthe,costhe,sinphi,cosphi;
            float r11, r12, r13, r21, r22, r23, r31, r32, r33;
			float a11, a33, a44, a66, a11a66, a13a44;

			vp2=vp0[k][i][j]*vp0[k][i][j];
            vs2=vs0[k][i][j]*vs0[k][i][j];
            ep=1+2*epsilon[k][i][j];
            de=1+2*delta[k][i][j];
            ga=1+2*gama[k][i][j];
            the=theta[k][i][j];
            phi=phai[k][i][j];

            costhe=cos(the);
            sinthe=sin(the);
            cosphi=cos(phi);
            sinphi=sin(phi);

            r11=costhe*cosphi;
            r12=-sinphi;
            r13=-sinthe*cosphi;
            r21=costhe*sinphi;
            r22=cosphi;
            r23=-sinthe*sinphi;
            r31=sinthe;
            r32=0.0;
            r33=costhe;

			a11=vp2*ep;
			a33=vp2;
			a44=vs2;
			a66=vs2*ga;
            vpn2=vp2*de;
            a11a66=a11-a66;
            a13a44=sqrt((vp2-vs2)*(vpn2-vs2));

          float pxy, pxz, pyz, qxy, qxz, qyz, rxy, rxz, ryz;
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
          float hpy, hqy, hry, hpx, hqx, hrx, hpz, hqz, hrz;
          // 2nd-order derivatives
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

		float px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2;
	    float qxy1, rxz1, pxy1, ryz1, pxz1, qyz1;

          // x-component
          px2 = r11*r11*hpx+r21*r21*hpy+r31*r31*hpz+2*r11*r21*pxy+2*r11*r31*pxz+2*r21*r31*pyz;
          py2 = r12*r12*hpx+r22*r22*hpy+r32*r32*hpz+2*r12*r22*pxy+2*r12*r32*pxz+2*r22*r32*pyz;
          pz2 = r13*r13*hpx+r23*r23*hpy+r33*r33*hpz+2*r13*r23*pxy+2*r13*r33*pxz+2*r23*r33*pyz;

          qxy1 = r11*r12*hqx+r21*r22*hqy+r31*r32*hqz+(r11*r22+r21*r12)*qxy+(r11*r32+r31*r12)*qxz+(r21*r32+r31*r22)*qyz;
          rxz1 = r11*r13*hrx+r21*r23*hry+r31*r33*hrz+(r11*r23+r21*r13)*rxy+(r11*r33+r31*r13)*rxz+(r21*r33+r31*r23)*ryz;

          p3[k][i][j]=2*p2[k][i][j] - p1[k][i][j] + dt2*( a11*px2 +  a66*py2 + a44*pz2 + a11a66*qxy1 + a13a44*rxz1);

          // y-component
          qx2 = r11*r11*hqx+r21*r21*hqy+r31*r31*hqz+2*r11*r21*qxy+2*r11*r31*qxz+2*r21*r31*qyz;
          qy2 = r12*r12*hqx+r22*r22*hqy+r32*r32*hqz+2*r12*r22*qxy+2*r12*r32*qxz+2*r22*r32*qyz;
          qz2 = r13*r13*hqx+r23*r23*hqy+r33*r33*hqz+2*r13*r23*qxy+2*r13*r33*qxz+2*r23*r33*qyz;

          pxy1 = r11*r12*hpx+r21*r22*hpy+r31*r32*hpz+(r11*r22+r21*r12)*pxy+(r11*r32+r31*r12)*pxz+(r21*r32+r31*r22)*pyz;
          ryz1 = r12*r13*hrx+r22*r23*hry+r32*r33*hrz+(r12*r23+r22*r13)*rxy+(r12*r33+r32*r13)*rxz+(r22*r33+r32*r23)*ryz;

          q3[k][i][j]=2*q2[k][i][j] - q1[k][i][j] + dt2*( a66*qx2 +  a11*qy2 + a44*qz2 + a11a66*pxy1 + a13a44*ryz1);

          // z-component
          rx2 = r11*r11*hrx+r21*r21*hry+r31*r31*hrz+2*r11*r21*rxy+2*r11*r31*rxz+2*r21*r31*ryz;
          ry2 = r12*r12*hrx+r22*r22*hry+r32*r32*hrz+2*r12*r22*rxy+2*r12*r32*rxz+2*r22*r32*ryz;
          rz2 = r13*r13*hrx+r23*r23*hry+r33*r33*hrz+2*r13*r23*rxy+2*r13*r33*rxz+2*r23*r33*ryz;

          pxz1 = r11*r13*hpx+r21*r23*hpy+r31*r33*hpz+(r11*r23+r21*r13)*pxy+(r11*r33+r31*r13)*pxz+(r21*r33+r31*r23)*pyz;
          qyz1 = r12*r13*hqx+r22*r23*hqy+r32*r33*hqz+(r12*r23+r22*r13)*qxy+(r12*r33+r32*r13)*qxz+(r22*r33+r32*r23)*qyz;

          r3[k][i][j]=2*r2[k][i][j] - r1[k][i][j] + dt2*( a44*rx2 +  a44*ry2 + a33*rz2 + a13a44*pxz1 + a13a44*qyz1);
	}
}

void fwpttielastic3dhomo(float dt2,float***p1,float***p2,float***p3,float***q1,float***q2,float***q3,
		             float***r1,float***r2,float***r3,float***px_tmp,float***qy_tmp,float***rz_tmp,
                     float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,
					 float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                     float dx, float dy, float dz, int nxpad, int nypad, int nzpad, 
                     float vp0,float vs0, float epsilon, float delta, float gama,
					 float theta, float phai)
/*< fwpttielastic3dhomo: forward-propagating using original elastic equation of displacement in 3D homogeneous TTI media>*/
{
    int   i,j,k,l;
    float r11, r12, r13, r21, r22, r23, r31, r32, r33;
	float qxy1, rxz1, pxy1, ryz1, pxz1, qyz1;
    float hpy, hqy, hry, hpx, hqx, hrx, hpz, hqz, hrz;
    float pxy, pxz, pyz, qxy, qxz, qyz, rxy, rxz, ryz;
		
	zero3float(px_tmp,nzpad,nxpad,nypad);	
	zero3float(qy_tmp,nzpad,nxpad,nypad);	
	zero3float(rz_tmp,nzpad,nxpad,nypad);	

#pragma omp parallel for private(i,j,k,l) \
	    schedule(dynamic) \
        shared(p2,q2,r2,px_tmp,qy_tmp,rz_tmp,coeff_1dx,coeff_1dy,coeff_1dz,nxpad,nypad,nzpad,dx,dy,dz)
	for(k=0;k<nypad;k++)
		for(i=0;i<nxpad;i++)
			for(j=0;j<nzpad;j++)
				for(l=-_mix;l<=_mix;l++)
				{
					if(i+l>=0&&i+l<nxpad)
						px_tmp[k][i][j]+=coeff_1dx[l+_mix]*p2[k][i+l][j]/2.0/dx;
					if(k+l>=0&&k+l<nypad)
						qy_tmp[k][i][j]+=coeff_1dy[l+_mix]*q2[k+l][i][j]/2.0/dy;
					if(j+l>=0&&j+l<nzpad)
						rz_tmp[k][i][j]+=coeff_1dz[l+_mix]*r2[k][i][j+l]/2.0/dz;
				}

#pragma omp parallel for private(i,j,k,l) \
	    schedule(dynamic) \
	    shared(p1,p2,p3,q1,q2,q3,r1,r2,r3,\
		px_tmp,qy_tmp,rz_tmp,\
		coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz,\
	    vp0,vs0,epsilon,delta,gama,theta,phai)
	for(k=0;k<nypad;k++)
	for(i=0;i<nxpad;i++)
		for(j=0;j<nzpad;j++)
		{
			float vp2,vs2,ep,de,ga,vpn2;
			float sinthe,costhe,sinphi,cosphi;
			float a11, a33, a44, a66, a11a66, a13a44;
			float px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2;
			vp2=vp0*vp0;
            vs2=vs0*vs0;
            ep=1+2*epsilon;
            de=1+2*delta;
            ga=1+2*gama;

            costhe=cos(theta);
            sinthe=sin(theta);
            cosphi=cos(phai);
            sinphi=sin(phai);

            r11=costhe*cosphi;
            r12=-sinphi;
            r13=-sinthe*cosphi;
            r21=costhe*sinphi;
            r22=cosphi;
            r23=-sinthe*sinphi;
            r31=sinthe;
            r32=0.0;
            r33=costhe;

			a11=vp2*ep;
			a33=vp2;
			a44=vs2;
			a66=vs2*ga;
            vpn2=vp2*de;
            a11a66=a11-a66;
            a13a44=sqrt((vp2-vs2)*(vpn2-vs2));

		  pxy=0;
		  pxz=0;
		  qxy=0;
		  qyz=0;
		  rxz=0;
		  ryz=0;
		  for(l=-_mix;l<=_mix;l++)
		  {
			  if(k+l>=0&&k+l<nypad)
			  {
				  pxy+=coeff_1dy[l+_mix]*px_tmp[k+l][i][j]/2.0/dy;
				  ryz+=coeff_1dz[l+_mix]*rz_tmp[k+l][i][j]/2.0/dy;
			  }
			  if(i+l>=0&&i+l<nxpad)
			  {
				  qxy+=coeff_1dy[l+_mix]*qy_tmp[k][i+l][j]/2.0/dx;
				  rxz+=coeff_1dz[l+_mix]*rz_tmp[k][i+l][j]/2.0/dx;
			  }
			  if(j+l>=0&&j+l<nzpad)
			  {
				  pxz+=coeff_1dz[l+_mix]*px_tmp[k][i][j+l]/2.0/dz;
				  qyz+=coeff_1dz[l+_mix]*qy_tmp[k][i][j+l]/2.0/dz;
			  }
		  }

          // 2nd-order derivatives
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

          // x-component
          px2 = r11*r11*hpx+r21*r21*hpy+r31*r31*hpz+2*r11*r21*pxy+2*r11*r31*pxz+2*r21*r31*pyz;
          py2 = r12*r12*hpx+r22*r22*hpy+r32*r32*hpz+2*r12*r22*pxy+2*r12*r32*pxz+2*r22*r32*pyz;
          pz2 = r13*r13*hpx+r23*r23*hpy+r33*r33*hpz+2*r13*r23*pxy+2*r13*r33*pxz+2*r23*r33*pyz;

          qxy1 = r11*r12*hqx+r21*r22*hqy+r31*r32*hqz+(r11*r22+r21*r12)*qxy+(r11*r32+r31*r12)*qxz+(r21*r32+r31*r22)*qyz;
          rxz1 = r11*r13*hrx+r21*r23*hry+r31*r33*hrz+(r11*r23+r21*r13)*rxy+(r11*r33+r31*r13)*rxz+(r21*r33+r31*r23)*ryz;

          p3[k][i][j]=2*p2[k][i][j] - p1[k][i][j] + dt2*( a11*px2 +  a66*py2 + a44*pz2 + a11a66*qxy1 + a13a44*rxz1);

          // y-component
          qx2 = r11*r11*hqx+r21*r21*hqy+r31*r31*hqz+2*r11*r21*qxy+2*r11*r31*qxz+2*r21*r31*qyz;
          qy2 = r12*r12*hqx+r22*r22*hqy+r32*r32*hqz+2*r12*r22*qxy+2*r12*r32*qxz+2*r22*r32*qyz;
          qz2 = r13*r13*hqx+r23*r23*hqy+r33*r33*hqz+2*r13*r23*qxy+2*r13*r33*qxz+2*r23*r33*qyz;

          pxy1 = r11*r12*hpx+r21*r22*hpy+r31*r32*hpz+(r11*r22+r21*r12)*pxy+(r11*r32+r31*r12)*pxz+(r21*r32+r31*r22)*pyz;
          ryz1 = r12*r13*hrx+r22*r23*hry+r32*r33*hrz+(r12*r23+r22*r13)*rxy+(r12*r33+r32*r13)*rxz+(r22*r33+r32*r23)*ryz;

          q3[k][i][j]=2*q2[k][i][j] - q1[k][i][j] + dt2*( a66*qx2 +  a11*qy2 + a44*qz2 + a11a66*pxy1 + a13a44*ryz1);

          // z-component
          rx2 = r11*r11*hrx+r21*r21*hry+r31*r31*hrz+2*r11*r21*rxy+2*r11*r31*rxz+2*r21*r31*ryz;
          ry2 = r12*r12*hrx+r22*r22*hry+r32*r32*hrz+2*r12*r22*rxy+2*r12*r32*rxz+2*r22*r32*ryz;
          rz2 = r13*r13*hrx+r23*r23*hry+r33*r33*hrz+2*r13*r23*rxy+2*r13*r33*rxz+2*r23*r33*ryz;

          pxz1 = r11*r13*hpx+r21*r23*hpy+r31*r33*hpz+(r11*r23+r21*r13)*pxy+(r11*r33+r31*r13)*pxz+(r21*r33+r31*r23)*pyz;
          qyz1 = r12*r13*hqx+r22*r23*hqy+r32*r33*hqz+(r12*r23+r22*r13)*qxy+(r12*r33+r32*r13)*qxz+(r22*r33+r32*r23)*qyz;

          r3[k][i][j]=2*r2[k][i][j] - r1[k][i][j] + dt2*( a44*rx2 +  a44*ry2 + a33*rz2 + a13a44*pxz1 + a13a44*qyz1);
	}
}
