/*************************************************************************
 * * Forward propagating using original elastic equation of displacement 
 *   in VTI media
 * 
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "_cjb.h"
#include "_fd.h"

#include "zero.h"

void fwpvtielastic(float dt2, float** p1,float** p2,float** p3, float** q1,float** q2,float** q3,
		   float* coeff_2dx,float* coeff_2dz, float* coeff_1dx,float* coeff_1dz,
		   float dx, float dz, int nx, int nz, int nxpad, int nzpad, 
		   float **vp0,float **vs0, float **epsilon,float **delta)
/*< fwpvtielastic: forward-propagating using original elastic equation of displacement in VTI media>*/
{
    int   i,j,l;

    int im, jm, lm;
		
    float **px_tmp;
    float **qx_tmp;

    float vp2;
    float vs2;
    float ep;
    float de;

    float vpx2;
    float vpn2;
    float coef;

    float px;
    float pz;
    float qx;
    float qz;
    float pxz;
    float qxz;

    px_tmp=sf_floatalloc2(nzpad,nxpad);
    qx_tmp=sf_floatalloc2(nzpad,nxpad);

    zero2float(px_tmp,nzpad,nxpad);	
    zero2float(qx_tmp,nzpad,nxpad);	

#ifdef _OPENMP
#pragma omp parallel for private(i,j,l,lm)	\
    schedule(dynamic)				\
    shared(p2,q2,px_tmp,qx_tmp,coeff_1dx,dx)
#endif
    for(i=_m;i<nx+_m;i++)
	for(j=_m;j<nz+_m;j++)
	{
	    for(l=-_m;l<=_m;l++)
	    {
		lm=l+_m;
		px_tmp[i][j]+=coeff_1dx[lm]*p2[i+l][j]/2.0/dx;
		qx_tmp[i][j]+=coeff_1dx[lm]*q2[i+l][j]/2.0/dx;
	    }
	}

#ifdef _OPENMP
#pragma omp parallel for private(i,j,l,im,jm,lm,vp2,vs2,ep,de,vpx2,vpn2,coef,px,pz,qx,qz,pxz,qxz) \
    schedule(dynamic)				\
    shared(p1,p2,p3,q1,q2,q3,px_tmp,qx_tmp,	\
	   coeff_1dx,coeff_1dz,			\
	   coeff_2dx,coeff_2dz,			\
	   vp0, vs0, epsilon, delta, dt2)
#endif
    for(i=_m;i<nx+_m;i++)
    {
	im=i-_m;
	for(j=_m;j<nz+_m;j++)
	{
	    jm=j-_m;

	    vp2=vp0[im][jm]*vp0[im][jm];
	    vs2=vs0[im][jm]*vs0[im][jm];
	    ep=1+2*epsilon[im][jm];
	    de=1+2*delta[im][jm];

	    vpx2=vp2*ep;
	    vpn2=vp2*de;
	    coef=sqrt((vp2-vs2)*(vpn2-vs2));

	    px=0;
	    pz=0;
	    qx=0;
	    qz=0;
	    pxz=0;
	    qxz=0;
	    for(l=-_m;l<=_m;l++)
	    {
		lm=l+_m;
		px+=coeff_2dx[lm]*p2[i+l][j];
		qx+=coeff_2dx[lm]*q2[i+l][j];
		pz+=coeff_2dz[lm]*p2[i][j+l];
		qz+=coeff_2dz[lm]*q2[i][j+l];
		pxz+=coeff_1dz[lm]*px_tmp[i][j+l]/2.0/dz;
		qxz+=coeff_1dz[lm]*qx_tmp[i][j+l]/2.0/dz;
	    }

	    p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*( vpx2*px + vs2*pz + coef*qxz);
	    q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*( vs2*qx + vp2*qz + coef*pxz);
	}
    }

    free(*px_tmp);	
    free(*qx_tmp);	
}

void fwpvtielastic3d(float dt2,float***p1,float***p2,float***p3,float***q1,float***q2,float***q3,
		     float***r1,float***r2,float***r3,float***px_tmp,float***qy_tmp,float***rz_tmp,
                     float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,
		     float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                     float dx, float dy, float dz, int nxpad, int nypad, int nzpad, 
                     float ***vp0,float ***vs0, float ***epsilon, float ***delta, float ***gama)
/*< fwpvtielastic3d: forward-propagating using original elastic equation of displacement in 3D VTI media>*/
{
    int   i,j,k,l;
    int il, jl, kl, lm;

    float vp2,vs2,ep,de,ga,vpn2;
    float a11, a33, a44, a66, a11a66, a13a44;
    float px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2;
    float pxy, pxz, qxy, qyz, rxz, ryz;

    zero3float(px_tmp,nzpad,nxpad,nypad);	
    zero3float(qy_tmp,nzpad,nxpad,nypad);	
    zero3float(rz_tmp,nzpad,nxpad,nypad);	

#ifdef _OPENMP
#pragma omp parallel for private(k,i,j,l)				\
    schedule(dynamic)							\
    shared(p2,q2,r2,px_tmp,qy_tmp,rz_tmp,coeff_1dx,coeff_1dy,coeff_1dz,nxpad,nypad,nzpad,dx,dy,dz)
#endif
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

#ifdef _OPENMP
#pragma omp parallel for private(k,i,j,l,vp2,vs2,ep,de,ga,vpn2,a11, a33, a44, a66, a11a66, a13a44,px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2, pxy, pxz, qxy, qyz, rxz, ryz, il, kl, jl, lm) \
    schedule(dynamic)						\
    shared(p1,p2,p3,q1,q2,q3,r1,r2,r3,px_tmp,qy_tmp,rz_tmp,	\
	   coeff_1dx,coeff_1dy,coeff_1dz,			\
	   coeff_2dx,coeff_2dy,coeff_2dz,			\
	   vp0, vs0, epsilon, delta, gama,dt2)
#endif
    for(k=0;k<nypad;k++)
	for(i=0;i<nxpad;i++)
	    for(j=0;j<nzpad;j++)
	    {

		vp2=vp0[k][i][j]*vp0[k][i][j];
		vs2=vs0[k][i][j]*vs0[k][i][j];
		ep=1+2*epsilon[k][i][j];
		de=1+2*delta[k][i][j];
		ga=1+2*gama[k][i][j];

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
		    kl=k+l;
		    if(kl>=0&&kl<nypad)
		    {
			pxy+=coeff_1dy[l+_mix]*px_tmp[kl][i][j]/2.0/dy;
			ryz+=coeff_1dz[l+_mix]*rz_tmp[kl][i][j]/2.0/dy;
		    }
		    il=i+l;
		    if(il>=0&&il<nxpad)
		    {
			qxy+=coeff_1dy[l+_mix]*qy_tmp[k][il][j]/2.0/dx;
			rxz+=coeff_1dz[l+_mix]*rz_tmp[k][il][j]/2.0/dx;
		    }
		    jl=j+l;
		    if(j+l>=0&&j+l<nzpad)
		    {
			pxz+=coeff_1dz[l+_mix]*px_tmp[k][i][jl]/2.0/dz;
			qyz+=coeff_1dz[l+_mix]*qy_tmp[k][i][jl]/2.0/dz;
		    }
		}

		/* 2nd-order derivatives */
		py2=0;
		qy2=0;
		ry2=0;
		px2=0;
		qx2=0;
		rx2=0;
		pz2=0;
		qz2=0;
		rz2=0;
		for(l=-_m;l<=_m;l++)
		{
		    lm=l+_m;
		    kl=k+l;
		    if(k+l>=0&&k+l<nypad)
		    {
			py2+=coeff_2dy[lm]*p2[kl][i][j];
			qy2+=coeff_2dy[lm]*q2[kl][i][j];
			ry2+=coeff_2dy[lm]*r2[kl][i][j];
		    }
		    il=i+l;
		    if(i+l>=0&&i+l<nxpad)
		    {
			px2+=coeff_2dx[lm]*p2[k][il][j];
			qx2+=coeff_2dx[lm]*q2[k][il][j];
			rx2+=coeff_2dx[lm]*r2[k][il][j];
		    }
		    jl=j+l;
		    if(j+l>=0&&j+l<nzpad)
		    {
			pz2+=coeff_2dz[lm]*p2[k][i][jl];
			qz2+=coeff_2dz[lm]*q2[k][i][jl];
			rz2+=coeff_2dz[lm]*r2[k][i][jl];
		    }
		}

		/* x-component */
		p3[k][i][j]=2*p2[k][i][j] - p1[k][i][j] + dt2*( a11*px2 +  a66*py2 + a44*pz2 + a11a66*qxy + a13a44*rxz);

		/* y-component */
		q3[k][i][j]=2*q2[k][i][j] - q1[k][i][j] + dt2*( a66*qx2 +  a11*qy2 + a44*qz2 + a11a66*pxy + a13a44*ryz);

		/* z-component */
		r3[k][i][j]=2*r2[k][i][j] - r1[k][i][j] + dt2*( a44*rx2 +  a44*ry2 + a33*rz2 + a13a44*pxz + a13a44*qyz);
	    }
}

void fwpvtielastic3dhomo(float dt2,float***p1,float***p2,float***p3,float***q1,float***q2,float***q3,
			 float***r1,float***r2,float***r3,float***px_tmp,float***qy_tmp,float***rz_tmp,
			 float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,
			 float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
			 float dx, float dy, float dz, int nx, int ny, int nz, 
			 float vp0,float vs0, float epsilon, float delta, float gama)
/*< fwpvtielastic3dhomo: forward-propagating using original elastic equation of displacement in 3D homogeneous VTI media>*/
{
    int   i,j,k,l;
    float vp2,vs2,ep,de,ga,vpn2;
    float a11, a33, a44, a66, a11a66, a13a44;
    float px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2;
    float pxy, pxz, qxy, qyz, rxz, ryz;

    int il, jl, kl, lm;

    zero3float(px_tmp,nz,nx,ny);	
    zero3float(qy_tmp,nz,nx,ny);	
    zero3float(rz_tmp,nz,nx,ny);	

#ifdef _OPENMP
#pragma omp parallel for private(k,i,j,l)				\
    schedule(dynamic)							\
    shared(p2,q2,r2,px_tmp,qy_tmp,rz_tmp,coeff_1dx,coeff_1dy,coeff_1dz,nx,ny,nz,dx,dy,dz)
#endif
    for(k=0;k<ny;k++)
	for(i=0;i<nx;i++)
	    for(j=0;j<nz;j++)
		for(l=-_mix;l<=_mix;l++)
		{
		    if(i+l>=0&&i+l<nx)
			px_tmp[k][i][j]+=coeff_1dx[l+_mix]*p2[k][i+l][j]/2.0/dx;
		    if(k+l>=0&&k+l<ny)
			qy_tmp[k][i][j]+=coeff_1dy[l+_mix]*q2[k+l][i][j]/2.0/dy;
		    if(j+l>=0&&j+l<nz)
			rz_tmp[k][i][j]+=coeff_1dz[l+_mix]*r2[k][i][j+l]/2.0/dz;
		}

#ifdef _OPENMP
#pragma omp parallel for private(k,i,j,l,vp2,vs2,ep,de,ga,vpn2,a11, a33, a44, a66, a11a66, a13a44, px2, py2, pz2, qx2, qy2, qz2, rx2, ry2, rz2, pxy, pxz, qxy, qyz, rxz, ryz, kl, il, jl, lm) \
    schedule(dynamic)						\
    shared(p1,p2,p3,q1,q2,q3,r1,r2,r3,px_tmp,qy_tmp,rz_tmp,	\
	   coeff_1dx,coeff_1dy,coeff_1dz,			\
	   coeff_2dx,coeff_2dy,coeff_2dz,			\
	   vp0, vs0, epsilon, delta, gama,dt2)
#endif
    for(k=0;k<ny;k++)
	for(i=0;i<nx;i++)
	    for(j=0;j<nz;j++)
	    {
		vp2=vp0*vp0;
		vs2=vs0*vs0;
		ep=1+2*epsilon;
		de=1+2*delta;
		ga=1+2*gama;

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
		    kl=k+l;
		    if(kl>=0&&kl<ny)
		    {
			pxy+=coeff_1dy[l+_mix]*px_tmp[kl][i][j]/2.0/dy;
			ryz+=coeff_1dz[l+_mix]*rz_tmp[kl][i][j]/2.0/dy;
		    }
		    il=i+l;
		    if(il>=0&&il<nx)
		    {
			qxy+=coeff_1dy[l+_mix]*qy_tmp[k][il][j]/2.0/dx;
			rxz+=coeff_1dz[l+_mix]*rz_tmp[k][il][j]/2.0/dx;
		    }
		    jl=j+l;
		    if(j+l>=0&&j+l<nz)
		    {
			pxz+=coeff_1dz[l+_mix]*px_tmp[k][i][jl]/2.0/dz;
			qyz+=coeff_1dz[l+_mix]*qy_tmp[k][i][jl]/2.0/dz;
		    }
		}

		/* 2nd-order derivatives */
		py2=0;
		qy2=0;
		ry2=0;
		px2=0;
		qx2=0;
		rx2=0;
		pz2=0;
		qz2=0;
		rz2=0;
		for(l=-_m;l<=_m;l++)
		{
		    lm=l+_m;
		    kl=k+l;
		    if(k+l>=0&&k+l<ny)
		    {
			py2+=coeff_2dy[lm]*p2[kl][i][j];
			qy2+=coeff_2dy[lm]*q2[kl][i][j];
			ry2+=coeff_2dy[lm]*r2[kl][i][j];
		    }
		    il=i+l;
		    if(i+l>=0&&i+l<nx)
		    {
			px2+=coeff_2dx[lm]*p2[k][il][j];
			qx2+=coeff_2dx[lm]*q2[k][il][j];
			rx2+=coeff_2dx[lm]*r2[k][il][j];
		    }
		    jl=j+l;
		    if(j+l>=0&&j+l<nz)
		    {
			pz2+=coeff_2dz[lm]*p2[k][i][jl];
			qz2+=coeff_2dz[lm]*q2[k][i][jl];
			rz2+=coeff_2dz[lm]*r2[k][i][jl];
		    }
		}

		/* x-component */
		p3[k][i][j]=2*p2[k][i][j] - p1[k][i][j] + dt2*( a11*px2 +  a66*py2 + a44*pz2 + a11a66*qxy + a13a44*rxz);

		/* y-component */
		q3[k][i][j]=2*q2[k][i][j] - q1[k][i][j] + dt2*( a66*qx2 +  a11*qy2 + a44*qz2 + a11a66*pxy + a13a44*ryz);

		/* z-component */
		r3[k][i][j]=2*r2[k][i][j] - r1[k][i][j] + dt2*( a44*rx2 +  a44*ry2 + a33*rz2 + a13a44*pxz + a13a44*qyz);
	    }
}
