/*
  Copyright (C) 2000 The Board of Trustees of Stanford University
  
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

#define EPS 1.0e-6

static void Lint3D(float xx,float yy,float zz,float val,float *model);

static int nx,ny,nz,nv;
static float ox,oy,oz,dx,dy,dz;

void RayTrace(int nx1,int ny1,int nz1,
	      float ox1, float oy1, float oz1,
	      float dx1, float dy1, float dz1,
	      float rx,float ry,float rz,
	      float sx,float sy,float sz,
	      float v0,float vgrad,
	      sf_file out)
/*< run ray tracing >*/
{
    int nxyz,i,ith,nth;
    float *output;
    double phi,theta,thetas,thetar;

    double z1,z2,x2, orx,ory;
    double beta,xr,rr,oth,dth;

    float xs,ys,zs,z0;

    nx=nx1; ny=ny1; nz=nz1; nv=nx*ny*nz;
    ox=ox1; oy=oy1; oz=oz1;
    dx=dx1; dy=dy1; dz=dz1;

    nth=500;
    nxyz=nx*ny*nz;
    output=sf_floatalloc(nxyz);
    for (i=0; i<nxyz; i++) output[i]=0.;

    if (vgrad > EPS) {

	x2=sqrt((sx-rx)*(sx-rx) + (sy-ry)*(sy-ry));
	phi=atan(fabs(sy-ry)/fabs(sx-rx));
	z0=v0/vgrad;
	if (sz<rz) {
	    z1=sz+z0; z2=rz+z0;
	    orx=sx;   ory=sy;
	} else {
	    z1=rz+z0; z2=sz+z0;
	    orx=rx;   ory=ry;
	}

	beta=0.5*(z1+z2)/x2;
	xr=0.5*x2-beta*(z1-z2);
	thetas=atan(z1/xr);
	thetar=atan(z2/(xr-x2));
	oth=thetas;
	dth=(thetar-thetas)/nth;
	rr=sqrt(xr*xr+z1*z1);

	for (ith=0; ith<nth; ith++) {
	    theta=oth+ith*dth;

	    xs=orx+(xr-rr*cos(theta))*cos(phi);
	    ys=ory+(xr-rr*cos(theta))*sin(phi);
      
	    zs=rr*sin(theta);
	    zs=zs-z0;

	    Lint3D(xs,ys,zs,1.,output);
	} 
    } else {
	for (ith=0; ith<nth; ith++) {
	    xs=sx+(rx-sx)*ith/(nth-1);
	    ys=sy+(ry-sy)*ith/(nth-1);
	    zs=sz+(rz-sz)*ith/(nth-1);
	    Lint3D(xs,ys,zs,1.,output);
	}
    }
    
    sf_floatwrite(output,nxyz,out);
}

static void Lint3D(float xx,float yy,float zz,float val,float *model)
{
  int ix,iy,iz,iv;
  float xc,yc,zc, fx,fy,fz;
  xc=(xx-ox)/dx; ix=(int) xc; fx=xc-ix;
  yc=(yy-oy)/dy; iy=(int) yc; fy=yc-iy;
  zc=(zz-oz)/dz; iz=(int) zc; fz=zc-iz;
  iv=(ix  )+(iy  )*nx+(iz  )*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(1-fx)*(1-fy)*(1-fz);
  iv=(ix+1)+(iy  )*nx+(iz  )*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(  fx)*(1-fy)*(1-fz);
  iv=(ix  )+(iy+1)*nx+(iz  )*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(1-fx)*(  fy)*(1-fz);
  iv=(ix+1)+(iy+1)*nx+(iz  )*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(  fx)*(  fy)*(1-fz);
  iv=(ix  )+(iy  )*nx+(iz+1)*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(1-fx)*(1-fy)*(  fz);
  iv=(ix+1)+(iy  )*nx+(iz+1)*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(  fx)*(1-fy)*(  fz);
  iv=(ix  )+(iy+1)*nx+(iz+1)*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(1-fx)*(  fy)*(  fz);
  iv=(ix+1)+(iy+1)*nx+(iz+1)*nx*ny;
  if ((iv>=0)&&(iv<nv)) model[iv]+=(  fx)*(  fy)*(  fz);
  return;
}

