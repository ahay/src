/*************************************************************************
 * * Forward propagating using elastic wave equation in ORT media
 * *
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
#include "_cjb.h"
#include "_fd.h"
#include "zero.h"

void fwportelastic(float dt2,float*** p1,float*** p2,float*** p3,float*** q1,float*** q2,float*** q3,
                   float*** r1,float*** r2,float*** r3,
                   float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                   float*** vp0,float ***vs0,float*** epsi_1,float*** del_1, float ***gama_1,float*** epsi_2,float*** del_2,
                   float ***gama_2,float ***del_3, float ***alpha, float ***the, float ***phi,
                   int nx, int ny, int nz, int nxpad, int nypad, int nzpad, float dx, float dy, float dz)
/*< fwportelastic: forward-propagating in ORT media with elastic wave equation>*/
{
    int   i,j,k,l;

#pragma omp parallel for private(i,j,k,l) \
	schedule(dynamic) \
	shared(p1,p2,p3, \
	q1,q2,q3, \
	r1,r2,r3,\
	coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz, \
	vp0,vs0,epsi_1,del_1,gama_1,epsi_2,del_2,gama_2,del_3)
    for(j=0;j<ny;j++)
    for(i=0;i<nx;i++)
    for(k=0;k<nz;k++){
    double px,py,pz,qx,qy,qz,rx,ry,rz;
    double vp2,vs2;
    double ep_1,de_1,ga_1,ep_2,de_2,ga_2,de_3;
    double vpx,vpy,vpz,vsz1,vsz2,vsz3,vpn1,vpn2,vpn3;
    double C23_44,C12_66,C13_55;
	vp2=vp0[j][i][k]*vp0[j][i][k];
	//tmp=(epsi_1[j][i][k]-del_1[j][i][k]);
	//vs2=tmp*vp2;
       	vs2=vs0[j][i][k]*vs0[j][i][k];//*0.5*0.5;                   
//	vs2=0.6*0.6*vp2;
	ep_1=1+2*epsi_1[j][i][k];
	de_1=1+2*del_1[j][i][k];
	ga_1=1+2*gama_1[j][i][k];
	ep_2=1+2*epsi_2[j][i][k];
	de_2=1+2*del_2[j][i][k];
	de_3=1+2*del_3[j][i][k];
	ga_2=1+2*gama_2[j][i][k];

	vpz=vp2;
	vpx=vp2*ep_2;
	vpy=vp2*ep_1;
	vpn1=vp2*de_1;
	vpn2=vp2*de_2;
	vpn3=vpx*de_3;//vpx*del_3
	vsz1=vs2;
	vsz2=vs2*ga_1/ga_2;
	vsz3=vs2*ga_1;
	C23_44=sqrt(vpz-vsz2)*sqrt(vpn1-vsz2);
	C12_66=sqrt(vpx-vsz3)*sqrt(vpn3-vsz3);
	C13_55=sqrt(vpz-vsz1)*sqrt(vpn2-vsz1);
	//deri calculation
	px=0;py=0;pz=0;
	qx=0;qy=0;qz=0;
	rx=0;ry=0;rz=0;
	for(l=-_m;l<=_m;l++)
	{
		if(i+l>=0&&i+l<nxpad)
		{
			px+=coeff_2dx[l+_m]*p2[j][i+l][k];
			qx+=coeff_2dx[l+_m]*q2[j][i+l][k];
			rx+=coeff_2dx[l+_m]*r2[j][i+l][k];
		}
		if(j+l>=0&&j+l<nypad)
		{
			py+=coeff_2dy[l+_m]*p2[j+l][i][k];
			qy+=coeff_2dy[l+_m]*q2[j+l][i][k];
			ry+=coeff_2dy[l+_m]*r2[j+l][i][k];
		}
               	if(k+l>=0&&k+l<nzpad)
		{
			pz+=coeff_2dz[l+_m]*p2[j][i][k+l];
			qz+=coeff_2dz[l+_m]*q2[j][i][k+l];
			rz+=coeff_2dz[l+_m]*r2[j][i][k+l];
		}
	}
	p3[j][i][k]=2*p2[j][i][k] - p1[j][i][k] + (float)(dt2*(vpx*px + vsz3*py + vsz1*pz + C12_66*qx + C13_55*rx));
	q3[j][i][k]=2*q2[j][i][k] - q1[j][i][k] + (float)(dt2*(C12_66*py + vsz3*qx + vpy*qy + vsz2*qz + C23_44*ry));
	r3[j][i][k]=2*r2[j][i][k] - r1[j][i][k] + (float)(dt2*(C13_55*pz + C23_44*qz + vsz1*rx + vsz2*ry + vpz*rz));
    }
}
