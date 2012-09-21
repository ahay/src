/*************************************************************************
* Forward propagating using pseudo-pure P-wave equation in VTI media
* (see, Kang and Cheng, 2011; Cheng et al., 2012).
*
*    Copyright: Tongji University (Jiubing Cheng)
*    2012.3.2
*************************************************************************/
#include <rsf.h>
#include "_cjb.h"
#include "_fd.h"
#include "alloc.h"

void fwpttipseudop(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, float dx,float dz,float dt,
               int nx,int nz,int nxpad, int nzpad, float** vp0,float **vs0,
               float** epsilon,float** delta, float **theta)
/*< fwpttipseudop: forward-propagating in TTI media with pseudo-pure P-wave equation>*/
{
        int i,j,k, im,jm,km;
        float **p_temp, **q_temp;
        float px,pz,qx,qz,vp2,vs2,vpx2,vpn2,ep,de,the,coef;
        float sinthe,costhe,cos2,sin2,sin2a,hxp,hxq,hzp,hzq,pxz,qxz;

        p_temp=alloc2float(nzpad,nxpad);
        q_temp=alloc2float(nzpad,nxpad);

        zero2float(p_temp,nzpad,nxpad);
        zero2float(q_temp,nzpad,nxpad);

        /* z-dreivative in mixed derivative when tilt angle nonzero */
        for(i=mix;i<nx+mix;i++)
                for(j=mix;j<nz+mix;j++)
                {
                        p_temp[i][j]=(p2[i][j+1]-p2[i][j-1])/2.0/dz;
                        q_temp[i][j]=(q2[i][j+1]-q2[i][j-1])/2.0/dz;
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

                        for(k=-m;k<=m;k++)
                        {
                                km=k+m;
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

                       p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*(vpx2*hxp + vs2*hzp + coef*hxq );

                       q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*(coef*hzp + vs2*hxq + vp2*hzq);
                }

           }/* i llop */
        free2float(p_temp);
        free2float(q_temp);
}
