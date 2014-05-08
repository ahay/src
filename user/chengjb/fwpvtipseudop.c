/*************************************************************************
* Forward propagating using pseudo-pure P-wave equation in VTI media
* (see, Kang and Cheng, 2011; Cheng et al., 2012).
*
*************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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


#include "_cjb.h"
#include "_fd.h"

void fwpvtipseudop(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, int nx,int nz,
               float** vp0,float **vs0,float** epsilon,float** delta)
/*< fwpvtipseudop: forward-propagating in VTI media with pseudo-pure P-wave equation>*/
{
        int i,j,k, im,jm,km;
        float px,pz,qx,qz,vp2,vs2,vpx2,vpn2,ep,de,coef;

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
                        coef=sqrt((vpn2-vs2)*(vp2-vs2));

                        px=0;
                        qx=0;
                        pz=0;
                        qz=0;

                        for(k=-_m;k<=_m;k++)
                        {
                                km=k+_m;
                                px+=coeff_x[km]*p2[i+k][j];
                                pz+=coeff_z[km]*p2[i][j+k];
                                qx+=coeff_x[km]*q2[i+k][j];
                                qz+=coeff_z[km]*q2[i][j+k];
                        }

                        p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*( vpx2*px + vs2*pz + coef*qx );

                        q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*( coef*pz + vs2*qx + vp2*qz);
                }

           }/* i llop */
}

void bwpvtipseudop(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, int nx,int nz,
               float** vp0,float **vs0,float** epsilon,float** delta)
/*< bwpvtipseudop: forward-propagating in VTI media with pseudo-pure P-wave equation>*/
{
        int i,j,k, im,jm,km;
        float px,pz,qx,qz,vp2,vs2,vpx2,vpn2,ep,de,coef;

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
                        coef=sqrt((vpn2-vs2)*(vp2-vs2));
                        
                        px=0;
                        qx=0;
                        pz=0;
                        qz=0;

                        for(k=-_m;k<=_m;k++)
                        {
                                km=k+_m;
                                px+=coeff_x[km]*p2[i+k][j];
                                pz+=coeff_z[km]*p2[i][j+k];
                                qx+=coeff_x[km]*q2[i+k][j];
                                qz+=coeff_z[km]*q2[i][j+k];
                        }

                        p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*( vpx2*px + vs2*pz + coef*qx );

                        q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*( coef*pz + vs2*qx + vp2*qz);
                }

           }/* i llop */
}
