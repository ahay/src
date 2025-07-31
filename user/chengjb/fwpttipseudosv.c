/*************************************************************************
* Forward propagating using pseudo-pure P-wave equation in TTI media
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


#include <rsf.h>
#include "_cjb.h"
#include "_fd.h"
#include "zero.h"

void fwpttipseudosv(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, float dx,float dz,
               int nx,int nz,int nxpad, int nzpad, float** vp0,float **vs0,
               float** epsilon,float** delta, float **theta)
/*< fwpttipseudosv: forward-propagating in TTI media with pseudo-pure SV-wave equation>*/
{
        int i,j,k, im,jm,km;
        float **p_temp, **q_temp;
        float px,pz,qx,qz,vp2,vs2,vpx2,vpn2,ep,de,the,coef;
        float sinthe,costhe,cos2,sin2,sin2a,hxp,hxq,hzp,hzq,pxz,qxz;

        p_temp=sf_floatalloc2(nzpad,nxpad);
        q_temp=sf_floatalloc2(nzpad,nxpad);

        zero2float(p_temp,nzpad,nxpad);
        zero2float(q_temp,nzpad,nxpad);

        /* z-dreivative in mixed derivative when tilt angle nonzero */
        for(i=_mix;i<nx+_mix;i++)
                for(j=_mix;j<nz+_mix;j++)
                {
                        p_temp[i][j]=(p2[i][j+1]-p2[i][j-1])/2.0/dz;
                        q_temp[i][j]=(q2[i][j+1]-q2[i][j-1])/2.0/dz;
                }

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
                        the=theta[im][jm];

                        the=theta[im][jm];
                        sinthe=sin(the);
                        costhe=cos(the);
                        cos2=costhe*costhe;
                        sin2=sinthe*sinthe;
                        sin2a=2*sinthe*costhe;

                        px=0;
                        qx=0;
                        pz=0;
                        qz=0;

                        vpx2=vp2*ep;
                        vpn2=vp2*de;
                        coef=sqrt((vpn2-vs2)*(vp2-vs2));
                        
                        //sf_warning("vp2=%f vs2=%f ep=%f de=%f",vp2,vs2,ep,de);

                        for(k=-_m;k<=_m;k++)
                        {
                                km=k+_m;
                                px+=coeff_x[km]*p2[i+k][j];
                                pz+=coeff_z[km]*p2[i][j+k];
                                qx+=coeff_x[km]*q2[i+k][j];
                                qz+=coeff_z[km]*q2[i][j+k];
                        }

                       /* x-dreivative in mixed derivative when tilt angle nonzero */
                       pxz=(p_temp[i+1][j]-p_temp[i-1][j])/2.0/dx;
                       qxz=(q_temp[i+1][j]-q_temp[i-1][j])/2.0/dx;

                       /* rotating according to the tilt angle */
                       hxp = cos2*px + sin2*pz + sin2a*pxz;
                       hxq = cos2*qx + sin2*qz + sin2a*qxz;

                       hzp = px + pz - hxp;
                       hzq = qx + qz - hxq;

                       p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*(vpx2*hxp + vs2*hzp - coef*hzq );

                       q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*(-coef*hxp + vs2*hxq + vp2*hzq);
                }

           }/* i llop */
        free(*p_temp);
        free(*q_temp);
}
