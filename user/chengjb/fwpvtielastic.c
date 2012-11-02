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

#include "_cjb.h"
#include "_fd.h"

#include "zero.h"

void fwpvtielastic(float dt2, float** p1,float** p2,float** p3, float** q1,float** q2,float** q3,
               float* coeff_2dx,float* coeff_2dz, float* coeff_1dx,float* coeff_1dz,
               float dx, float dz, int nx, int nz, int nxpad, int nzpad, 
               float **vp0,float **vs0, float **epsilon,float **delta)
/*< fwpvtielastic: forward-propagating using original elastic equation of displacement in VTI media>*/
{
        int   i,j,l, lm, im, jm;
        float px,pz,pxz, qxz,qx,qz;
        float vp2,vs2,ep,de,vpx2,vpn2,coef;
		
	float **px_tmp=sf_floatalloc2(nzpad,nxpad);
	float **qx_tmp=sf_floatalloc2(nzpad,nxpad);

	zero2float(px_tmp,nzpad,nxpad);	
	zero2float(qx_tmp,nzpad,nxpad);	

        for(i=m;i<nx+m;i++)
	for(j=m;j<nz+m;j++)
	{
		for(l=-m;l<=m;l++)
		{
                        lm=l+m;
			px_tmp[i][j]+=coeff_1dx[lm]*p2[i+l][j]/2.0/dx;
			qx_tmp[i][j]+=coeff_1dx[lm]*q2[i+l][j]/2.0/dx;
		}
	}

        for(i=m;i<nx+m;i++)
        {
           im=i-m;
	   for(j=m;j<nz+m;j++)
	   {
               jm=j-m;

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
		for(l=-m;l<=m;l++)
		{
                     lm=l+m;
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
