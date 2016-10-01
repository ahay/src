/*************************************************************************
* Forward propagating using pseudo-pure SV-wave equation in VTI media
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

void fwpvtipseudosv(float dt2,float** p1,float** p2,float** p3,float** q1,float** q2,float** q3,
               float* coeff_x,float* coeff_z, int nx,int nz,
               float** vp0,float **vs0,float** epsilon,float** delta)
/*< fwpvtipseudosv: forward-propagating in VTI media with pseudo-pure SV-wave equation>*/
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
                        
                        //sf_warning("vp2=%f vs2=%f ep=%f de=%f",vp2,vs2,ep,de);

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

                        p3[i][j]=2*p2[i][j] - p1[i][j] + dt2*( vpx2*px + vs2*pz - coef*qz );

                        q3[i][j]=2*q2[i][j] - q1[i][j] + dt2*(-coef*px + vs2*qx + vp2*qz);
                }

           }/* i llop */
}

/******************************for the 3D VTI media****************************/
void fwpvti3dpseudosv(float dt,float***p1,float***p2,float***p3, 
		float***r1,float***r2,float***r3, 
		float*coeff_2dx,float*coeff_2dy,float*coeff_2dz,float*coeff_1dx,float*coeff_1dy,float*coeff_1dz,
		float ***vp0, float ***vs0,float ***epsi,float ***del, 
		int nx, int ny, int nz, float dx, float dy, float dz)
/*< fwpvti3dpseudosv: forward-propagating in VTI media with pseudo-pure SV-wave equation>*/
{

	int i,j,k,l, m=_m;
#pragma omp parallel for private(i,j,k,l) \
	schedule(dynamic) \
	shared(p1,p2,p3, \
	r1,r2,r3,\
	coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz, \
	vp0,vs0,epsi,del)

	for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
				float px,py,pz,rx,ry,rz;
				float vp2,vs2,dt2;
				float vpx,vpn;
				float ep,de;
				float C23_44;

				dt2=dt*dt;
				vp2=vp0[i][j][k]*vp0[i][j][k];
				vs2=vs0[i][j][k]*vs0[i][j][k];
				ep=1+2*epsi[i][j][k];
				de=1+2*del[i][j][k];
				vpx=vp2*ep;
				vpn=vp2*de;
				C23_44=sqrt(vp2-vs2)*sqrt(vpn-vs2);
				
				//deri calculation
				px=0;py=0;pz=0;
				rx=0;ry=0;rz=0;
				/* Note the 3-D qSV-wave's pseudo-pure-mode wave equation has summed
				 * x- and y- component */
				for(l=-m;l<=m;l++)
				{
					if(i+l>=0&&i+l<nx)
					{
						px+=coeff_2dx[l+m]*p2[i+l][j][k];
						rx+=coeff_2dx[l+m]*r2[i+l][j][k];
					}
					if(j+l>=0&&j+l<ny)
					{

						py+=coeff_2dy[l+m]*p2[i][j+l][k];
						ry+=coeff_2dy[l+m]*r2[i][j+l][k];
					}
                    if(k+l>=0&&k+l<nz)
					{
						pz+=coeff_2dz[l+m]*p2[i][j][k+l];
						rz+=coeff_2dz[l+m]*r2[i][j][k+l];
					}
				}
				p3[i][j][k]=2*p2[i][j][k] - p1[i][j][k] + dt2*(vpx*(px+py) + vs2*pz - C23_44*rz);
				r3[i][j][k]=2*r2[i][j][k] - r1[i][j][k] + dt2*(-C23_44*(px+py) + vs2*(rx+ry) + vp2*rz);
			}
}

/******************************for the 3D VTI media****************************/
void fwpvti3dpseudosvhomo(float dt,float***p1,float***p2,float***p3, 
		float***r1,float***r2,float***r3, 
		float*coeff_2dx,float*coeff_2dy,float*coeff_2dz,float*coeff_1dx,float*coeff_1dy,float*coeff_1dz,
		float vp0, float vs0,float epsi,float del, 
		int nx, int ny, int nz, float dx, float dy, float dz)
/*< fwpvtipseudosvhomo: forward-propagating in VTI media with pseudo-pure SV-wave equation>*/
{

	int i,j,k,l, m=_m;
#pragma omp parallel for private(i,j,k,l) \
	schedule(dynamic) \
	shared(p1,p2,p3, \
	r1,r2,r3,\
	coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz, \
	vp0,vs0,epsi,del)

	for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
				float px,py,pz,rx,ry,rz;
				float vp2,vs2,dt2;
				float vpx,vpn;
				float ep,de;
				float C23_44;

				dt2=dt*dt;
				vp2=vp0*vp0;
				vs2=vs0*vs0;
				ep=1+2*epsi;
				de=1+2*del;
				vpx=vp2*ep;
				vpn=vp2*de;
				C23_44=sqrt(vp2-vs2)*sqrt(vpn-vs2);
				//C12_66=sqrt(vpx-vs2)*sqrt(vpn-vs2);
				
				//deri calculation
				px=0;py=0;pz=0;
				rx=0;ry=0;rz=0;
				/* Note the 3-D qSV-wave's pseudo-pure-mode wave equation has summed
				 * x- and y- component */
				for(l=-m;l<=m;l++)
				{
					if(i+l>=0&&i+l<nx)
					{
						px+=coeff_2dx[l+m]*p2[i+l][j][k];
						rx+=coeff_2dx[l+m]*r2[i+l][j][k];
					}
					if(j+l>=0&&j+l<ny)
					{

						py+=coeff_2dy[l+m]*p2[i][j+l][k];
						ry+=coeff_2dy[l+m]*r2[i][j+l][k];
					}
                    if(k+l>=0&&k+l<nz)
					{
						pz+=coeff_2dz[l+m]*p2[i][j][k+l];
						rz+=coeff_2dz[l+m]*r2[i][j][k+l];
					}
				}
				p3[i][j][k]=2*p2[i][j][k] - p1[i][j][k] + dt2*(vpx*(px+py) + vs2*pz - C23_44*rz);
				r3[i][j][k]=2*r2[i][j][k] - r1[i][j][k] + dt2*(-C23_44*(px+py) + vs2*(rx+ry) + vp2*rz);
			}
}
