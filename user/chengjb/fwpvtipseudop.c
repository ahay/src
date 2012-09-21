/*************************************************************************
* Forward propagating using pseudo-pure P-wave equation in VTI media
* (see, Kang and Cheng, 2011; Cheng et al., 2012).
*
*    Copyright: Tongji University (Jiubing Cheng)
*    2012.3.2
*************************************************************************/
#include "_cjb.h"
#include "_fd.h"

void fwpvtipseudop(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, float dx,float dz,float dt,
               int nx,int nz,int nxpad, int nzpad, float** vp0,float **vs0,float** epsilon,float** delta)
/*< fwpvtipseudop: forward-propagating in VTI media with pseudo-pure P-wave equation>*/
{
        int i,j,k, im,jm,km;
        float px,pz,qx,qz,vp2,vs2,vpx2,vpn2,ep,de,coef;

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
                        coef=sqrt((vpn2-vs2)*(vp2-vs2));
                        
                        //sf_warning("vp2=%f vs2=%f ep=%f de=%f",vp2,vs2,ep,de);

                        px=0;
                        qx=0;
                        pz=0;
                        qz=0;

                        for(k=-m;k<=m;k++)
                        {
                                km=k+m;
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
