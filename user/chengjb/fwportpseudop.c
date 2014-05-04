/*************************************************************************
 * * Forward propagating using pseudo-pure P-wave equation in ORT media
 * * (see, Cheng and Kang, 2012; Kang and Cheng, 2011; Wang et al., 2012).
 * *
 * *************************************************************************/
/*
  Copyright (C) 2012 Tongji University, Shanghai, China 
  Authors: Jiubing Cheng, Tengfei Wang and Wei Kang
     
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

void fwportpseudop1(float dt2,float*** p1,float*** p2,float*** p3,float*** q1,float*** q2,float*** q3,
		    float*** r1,float*** r2,float*** r3,
		    float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
		    float*** vp0,float ***vs0,float*** epsi_1,float*** del_1, float ***gama_1,float*** epsi_2,float*** del_2,
		    float ***gama_2,float ***del_3, float ***alpha, float ***the, float ***phi,
		    int nx, int ny, int nz, int nxpad, int nypad, int nzpad, float dx, float dy, float dz)
/*< fwportpseudop1: forward-propagating in ORT media with pseudo-pure P-wave equation>*/
{
    int   i,j,k,l,lmix,lm,il,jl,kl,ii,jj,kk;
    float px,py,pz,qx,qy,qz,rx,ry,rz;
    float hpx,hpy,hpz,hqx,hqy,hqz,hrx,hry,hrz;
    float hpxz,hpxy,hpyz,hqxz,hqxy,hqyz,hrxz,hrxy,hryz;

    float vp2,vs2;
    float vpx,vpy,vpz,vsz1,vsz2,vsz3,vpn1,vpn2,vpn3;
    float cos_alpha,sin_alpha,cos_theta,sin_theta,cos_phi,sin_phi;
    float ep_1,de_1,gam_1,ep_2,de_2,gam_2,de_3;
        
    float RT00,RT01,RT02,RT10,RT11,RT12,RT20,RT21,RT22;

    float ***px_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float ***py_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float ***qx_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float ***qy_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float ***rx_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float ***ry_tmp=sf_floatalloc3(nzpad,nxpad,nypad);

    zero3float(px_tmp,nzpad,nxpad,nypad);
    zero3float(py_tmp,nzpad,nxpad,nypad);
    zero3float(qx_tmp,nzpad,nxpad,nypad);
    zero3float(qy_tmp,nzpad,nxpad,nypad);
    zero3float(rx_tmp,nzpad,nxpad,nypad);
    zero3float(ry_tmp,nzpad,nxpad,nypad);

    float r2dx,r2dy,r2dz;
    r2dx=0.5/dx;
    r2dy=0.5/dy;
    r2dz=0.5/dz;

    //sf_warning("r2dx=%f r2dy=%f r2dz=%f",r2dx,r2dy,r2dz);

    //prepare for mixed deri calculation
    for(i=_mix;i<ny+_mix;i++)
	for(j=_mix;j<nx+_mix;j++)
	    for(k=_mix;k<nz+_mix;k++)
	    {
		for(l=-_mix;l<=_mix;l++)
		{
		    lmix=l+_mix;
		    il=i+l;
		    jl=j+l;
		    px_tmp[i][j][k]+=coeff_1dx[lmix]*p2[il][j][k]*r2dx;
		    qx_tmp[i][j][k]+=coeff_1dx[lmix]*q2[il][j][k]*r2dx;
		    rx_tmp[i][j][k]+=coeff_1dx[lmix]*r2[il][j][k]*r2dx;

		    py_tmp[i][j][k]+=coeff_1dy[lmix]*p2[i][jl][k]*r2dy;
		    qy_tmp[i][j][k]+=coeff_1dy[lmix]*q2[i][jl][k]*r2dy;
		    ry_tmp[i][j][k]+=coeff_1dy[lmix]*r2[i][jl][k]*r2dy;
		}
	    }

    for(i=0;i<nypad;i++)
    {
	ii=i-_m;
	if(ii<0)   ii=0;
	if(ii>=ny) ii=ny-1;
	for(j=0;j<nxpad;j++)
	{
	    jj=j-_m;
	    if(jj<0)   jj=0;
	    if(jj>=nx) jj=nx-1;
	    for(k=0;k<nzpad;k++)
	    {
                kk=k-_m;
                if(kk<0)   kk=0;
                if(kk>=nz) kk=nz-1;
		vp2=vp0[ii][jj][kk]*vp0[ii][jj][kk];
                //float tmp;
                //tmp=(epsi_1[ii][jj][kk]-del_1[ii][jj][kk]);
                //vs2=tmp*vp2;
                //vs2=0.6*0.6*vp2;
                vs2=vs0[ii][jj][kk]*vs0[ii][jj][kk];//*0.5*0.5;                   
		ep_1=1+2*epsi_1[ii][jj][kk];
		ep_2=1+2*epsi_2[ii][jj][kk];
		de_1=1+2*del_1[ii][jj][kk];
		de_2=1+2*del_2[ii][jj][kk];
		de_3=1+2*del_3[ii][jj][kk];
		gam_1=1+2*gama_1[ii][jj][kk];
		gam_2=1+2*gama_2[ii][jj][kk];

	        cos_alpha=cos(alpha[ii][jj][kk]);
        	sin_alpha=sin(alpha[ii][jj][kk]);
        	cos_theta=cos(the[ii][jj][kk]);
        	sin_theta=sin(the[ii][jj][kk]);
        	cos_phi=cos(phi[ii][jj][kk]);
        	sin_phi=sin(phi[ii][jj][kk]);

	        RT00=cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi;
	        RT01=sin_alpha*cos_theta*cos_phi+cos_alpha*sin_phi;
	        RT02=sin_theta*cos_phi;
	        RT10=-cos_alpha*cos_theta*sin_phi-sin_alpha*cos_phi;
	        RT11=-sin_alpha*cos_theta*sin_phi+cos_alpha*cos_phi;
	        RT12=-sin_theta*sin_phi;
	        RT20=-cos_alpha*sin_theta;
	        RT21=-sin_alpha*sin_theta;
	        RT22=cos_theta;

                //sf_warning("vp0=%f vs0=%f ",vp0[ii][jj][kk],vs0[ii][jj][kk]);
                //sf_warning("ep1=%f ep2=%f ",epsi_1[ii][jj][kk],epsi_2[ii][jj][kk]);
                //sf_warning("de1=%f de2=%f de3=%f ",del_1[ii][jj][kk],del_2[ii][jj][kk],del_3[ii][jj][kk]);
                //sf_warning("ga1=%f ga2=%f ",gama_1[ii][jj][kk],gama_2[ii][jj][kk]);
                //sf_warning("alpha=%f the=%f phi=%f ",alpha[ii][jj][kk],the[ii][jj][kk],phi[ii][jj][kk]);

		vpz=vp2;
		vpx=vp2*ep_2;
		vpy=vp2*ep_1;
		vpn1=vp2*de_1;
		vpn2=vp2*de_2;
		vpn3=vpx*de_3;//vpx*del_3
		vsz1=vs2;
		vsz2=vs2*gam_1/gam_2;
		vsz3=vs2*gam_1;

		//mixed deri calculation
		hpxy=0;hpxz=0;hpyz=0;
		hqxy=0;hqxz=0;hqyz=0;
		hrxy=0;hrxz=0;hryz=0;
		for(l=-_mix;l<=_mix;l++)
                {
		    lmix=l+_mix;
		    jl=j+l;
		    kl=k+l;
		    if(jl>=0&&jl<nypad)
		    {
                        hpxy+=coeff_1dy[lmix]*px_tmp[i][jl][k]*r2dy;
                        hqxy+=coeff_1dy[lmix]*qx_tmp[i][jl][k]*r2dy;
                        hrxy+=coeff_1dy[lmix]*rx_tmp[i][jl][k]*r2dy;
		    }
		    if(kl>=0&&kl<nzpad)
		    {
                        hpxz+=coeff_1dz[lmix]*px_tmp[i][j][kl]*r2dz;
                        hqxz+=coeff_1dz[lmix]*qx_tmp[i][j][kl]*r2dz;
                        hrxz+=coeff_1dz[lmix]*rx_tmp[i][j][kl]*r2dz;
                        hpyz+=coeff_1dz[lmix]*py_tmp[i][j][kl]*r2dz;
                        hqyz+=coeff_1dz[lmix]*qy_tmp[i][j][kl]*r2dz;
                        hryz+=coeff_1dz[lmix]*ry_tmp[i][j][kl]*r2dz;
                    }
                }// l loop

		//deri calculation
		hpx=0;hpy=0;hpz=0;
		hqx=0;hqy=0;hqz=0;
		hrx=0;hry=0;hrz=0;
		for(l=-_m;l<=_m;l++)
		{
		    lm=l+_m;
		    il=i+l;
		    jl=j+l;
		    kl=k+l;
		    if(il>=0&&il<nxpad)
		    {
			hpx+=coeff_2dx[lm]*p2[il][j][k];
			hqx+=coeff_2dx[lm]*q2[il][j][k];
			hrx+=coeff_2dx[lm]*r2[il][j][k];
		    }
		    if(jl>=0&&jl<nypad)
		    {
			hpy+=coeff_2dy[lm]*p2[i][jl][k];
			hqy+=coeff_2dy[lm]*q2[i][jl][k];
			hry+=coeff_2dy[lm]*r2[i][jl][k];
		    }
		    if(kl>=0&&kl<nzpad)
		    {
			hpz+=coeff_2dz[lm]*p2[i][j][kl];
			hqz+=coeff_2dz[lm]*q2[i][j][kl];
			hrz+=coeff_2dz[lm]*r2[i][j][kl];
		    }
		}// l loop
		// Tilted deri calculation
		px=0;py=0;pz=0;
		qx=0;qy=0;qz=0;
		rx=0;ry=0;rz=0;

		px=RT00*RT00*hpx+RT01*RT01*hpy+RT02*RT02*hpz+2*RT00*RT01*hpxy+2*RT00*RT02*hpxz+2*RT01*RT02*hpyz;
		py=RT10*RT10*hpx+RT11*RT11*hpy+RT12*RT12*hpz+2*RT10*RT11*hpxy+2*RT10*RT12*hpxz+2*RT11*RT12*hpyz;
		pz=RT20*RT20*hpx+RT21*RT21*hpy+RT22*RT22*hpz+2*RT20*RT21*hpxy+2*RT20*RT22*hpxz+2*RT21*RT22*hpyz;

		qx=RT00*RT00*hqx+RT01*RT01*hqy+RT02*RT02*hqz+2*RT00*RT01*hqxy+2*RT00*RT02*hqxz+2*RT01*RT02*hqyz;
		qy=RT10*RT10*hqx+RT11*RT11*hqy+RT12*RT12*hqz+2*RT10*RT11*hqxy+2*RT10*RT12*hqxz+2*RT11*RT12*hqyz;
		qz=RT20*RT20*hqx+RT21*RT21*hqy+RT22*RT22*hqz+2*RT20*RT21*hqxy+2*RT20*RT22*hqxz+2*RT21*RT22*hqyz;

		rx=RT00*RT00*hrx+RT01*RT01*hry+RT02*RT02*hrz+2*RT00*RT01*hrxy+2*RT00*RT02*hrxz+2*RT01*RT02*hryz;
		ry=RT10*RT10*hrx+RT11*RT11*hry+RT12*RT12*hrz+2*RT10*RT11*hrxy+2*RT10*RT12*hrxz+2*RT11*RT12*hryz;
		rz=RT20*RT20*hrx+RT21*RT21*hry+RT22*RT22*hrz+2*RT20*RT21*hrxy+2*RT20*RT22*hrxz+2*RT21*RT22*hryz;

		p3[i][j][k]=2*p2[i][j][k] - p1[i][j][k] + dt2*( vpx*px + vsz3*py + vsz1*pz + sqrt((vpx-vsz3)*(vpn3-vsz3))*qx + sqrt((vpz-vsz1)*(vpn2-vsz1))*rx);
		q3[i][j][k]=2*q2[i][j][k] - q1[i][j][k] + dt2*(sqrt((vpx-vsz3)*(vpn3-vsz3))*py + vsz3*qx + vpy*qy + vsz2*qz + sqrt((vpz-vsz2)*(vpn1-vsz2))*ry);
		r3[i][j][k]=2*r2[i][j][k] - r1[i][j][k] + dt2*(sqrt((vpz-vsz1)*(vpn2-vsz1))*pz+sqrt((vpz-vsz2)*(vpn1-vsz2))*qz+ vsz1*rx+vsz2*ry + vpz*rz);
/*
  if(i==m+ny/2&&j==m+nx/2&&k==m+nz/2)
  {
  sf_warning("p3= %f ",p3[i][j][k]);
  sf_warning("q3= %f ",q3[i][j][k]);
  sf_warning("r3= %f ",r3[i][j][k]);
  }
*/
	    }// k llop
	}// j llop
    }// i loop
    free(**px_tmp);
    free(**py_tmp);
    free(**qx_tmp);
    free(**qy_tmp);
    free(**rx_tmp);
    free(**ry_tmp);
}

void fwportpseudop(float dt2,float*** p1,float*** p2,float*** p3,float*** q1,float*** q2,float*** q3,
                   float*** r1,float*** r2,float*** r3,
                   float* coeff_2dx,float* coeff_2dy,float* coeff_2dz,float* coeff_1dx,float* coeff_1dy,float* coeff_1dz,
                   float*** vp0,float ***vs0,float*** epsi_1,float*** del_1, float ***gama_1,float*** epsi_2,float*** del_2,
                   float ***gama_2,float ***del_3, float ***alpha, float ***the, float ***phi,
                   int nx, int ny, int nz, int nxpad, int nypad, int nzpad, float dx, float dy, float dz)
/*< fwportpseudop: forward-propagating in ORT media with pseudo-pure P-wave equation>*/
{
    int   i,j,k,l;
    double px,py,pz,qx,qy,qz,rx,ry,rz;

    double vp2,vs2;
    double ep_1,de_1,ga_1,ep_2,de_2,ga_2,de_3;
    double vpx,vpy,vpz,vsz1,vsz2,vsz3,vpn1,vpn2,vpn3;
        
    double C23_44,C12_66,C13_55;

    for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	    for(k=0;k<nz;k++){
		vp2=vp0[j][i][k]*vp0[j][i][k];
		//tmp=(epsi_1[j][i][k]-del_1[j][i][k]);
		//vs2=tmp*vp2;
		//vs2=vs0[j][i][k]*vs0[j][i][k];//*0.5*0.5;                   
		vs2=0.6*0.6*vp2;
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

/******************************for the hti media****************************/
void fwportpseudophomo(float dt,float***p1,float***p2,float***p3, 
		       float***q1,float***q2,float***q3, 
		       float***r1,float***r2,float***r3, 
		       float*coeff_2dx,float*coeff_2dy,float*coeff_2dz,float*coeff_1dx,float*coeff_1dy,float*coeff_1dz,
		       float vp0, float vs0,float epsi1,float del1, float gama1,float epsi2,float del2, 
		       float gama2,float del3, int nx, int ny, int nz, float dx, float dy, float dz)
/*< fwportpseudophomo: forward-propagating in ORT media with pseudo-pure P-wave equation>*/
{
    float px,py,pz,qx,qy,qz,rx,ry,rz;
    float vp2,vs2,dt2;
    float vpx,vpy,vpz,vsz1,vsz2,vsz3,vpn1,vpn2,vpn3;
    float ep_1,de_1,gam_1,ep_2,de_2,gam_2,de_3;
    float C23_44,C12_66,C13_55;

    int i,j,k,l, m=_m;
#pragma omp parallel for private(i,j,k,l)				\
    schedule(dynamic)							\
    shared(p1,p2,p3,							\
	   q1,q2,q3,							\
	   r1,r2,r3,							\
	   coeff_1dx,coeff_1dy,coeff_1dz,coeff_2dx,coeff_2dy,coeff_2dz, \
	   vp0,vs0,epsi1,del1,gama1,epsi2,del2,gama2,del3)

    for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {

		dt2=dt*dt;
		vp2=vp0*vp0;
		vs2=vs0*vs0;
/*		float tmp;
		tmp=(epsi1-del1); */
		ep_1=1+2*epsi1;
		de_1=1+2*del1;
		gam_1=1+2*gama1;
		ep_2=1+2*epsi2;
		de_2=1+2*del2;
		de_3=1+2*del3;
		gam_2=1+2*gama2;

		vpz=vp2;
		vpx=vp2*ep_2;
		vpy=vp2*ep_1;
		vpn1=vp2*de_1;
		vpn2=vp2*de_2;
		vpn3=vpx*de_3;//vpx*del_3
		vsz1=vs2;
		vsz2=vs2*gam_1/gam_2;
		vsz3=vs2*gam_1;
		C23_44=sqrt(vpz-vsz2)*sqrt(vpn1-vsz2);
		C12_66=sqrt(vpx-vsz3)*sqrt(vpn3-vsz3);
		C13_55=sqrt(vpz-vsz1)*sqrt(vpn2-vsz1);
		//deri calculation
		px=0;py=0;pz=0;
		qx=0;qy=0;qz=0;
		rx=0;ry=0;rz=0;
		for(l=-m;l<=m;l++)
		{
		    if(i+l>=0&&i+l<nx)
		    {
			px+=coeff_2dx[l+m]*p2[i+l][j][k];
			qx+=coeff_2dx[l+m]*q2[i+l][j][k];
			rx+=coeff_2dx[l+m]*r2[i+l][j][k];
		    }
		    if(j+l>=0&&j+l<ny)
		    {

			py+=coeff_2dy[l+m]*p2[i][j+l][k];
			qy+=coeff_2dy[l+m]*q2[i][j+l][k];
			ry+=coeff_2dy[l+m]*r2[i][j+l][k];
		    }
                    if(k+l>=0&&k+l<nz)
		    {
			pz+=coeff_2dz[l+m]*p2[i][j][k+l];
			qz+=coeff_2dz[l+m]*q2[i][j][k+l];
			rz+=coeff_2dz[l+m]*r2[i][j][k+l];
		    }
		}
		p3[i][j][k]=2*p2[i][j][k] - p1[i][j][k] + dt2*(vpx*px + vsz3*py + vsz1*pz + C12_66*qx + C13_55*rx);
		q3[i][j][k]=2*q2[i][j][k] - q1[i][j][k] + dt2*(C12_66*py + vsz3*qx + vpy*qy + vsz2*qz + C23_44*ry);
		r3[i][j][k]=2*r2[i][j][k] - r1[i][j][k] + dt2*(C13_55*pz + C23_44*qz + vsz1*rx + vsz2*ry + vpz*rz);
	    }
}

void fwportpseudop2(float dt2,float*** p1,float*** p2,float*** p3,float*** q1,float*** q2,float*** q3, float***r1,float***r2,float***r3,
		    float*coeff_2dx,float*coeff_2dy,float*coeff_2dz,float*coeff_1dx,float*coeff_1dy,float*coeff_1dz,
		    float***vp0,float ***vs0,float***epsi_1,float***del_1, float ***gama_1,float***epsi_2,float***del_2, 
		    float ***gama_2,float ***delta_3,
                    int nx, int ny, int nz, int nxpad, int nypad, int nzpad)
{

    int i,j,k,l;

    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++)
	    for(k=0;k<nz;k++)
	    {
                //sf_warning("i%d j=%d k=%d",i,j,k);
		float px,py,pz,qx,qy,qz,rx,ry,rz;
		float vp2,vs2;
		float vpx,vpy,vpz,vsz1,vsz2,vsz3,vpn1,vpn2,vpn3;
		float ep_1,de_1,gam_1,ep_2,de_2,gam_2,de_3;
		float C23_44,C12_66,C13_55;

		vp2=vp0[i][j][k]*vp0[i][j][k];
		float tmp;
		tmp=(epsi_1[i][j][k]-del_1[i][j][k]);
		vs2=tmp*vp2;
		vs2=0.6*0.6*vp2;
                vs2=vs0[i][j][k]*vs0[i][j][k];//*0.5*0.5;                   
		ep_1=1+2*epsi_1[i][j][k];
		de_1=1+2*del_1[i][j][k];
		gam_1=1+2*gama_1[i][j][k];
		ep_2=1+2*epsi_2[i][j][k];
		de_2=1+2*del_2[i][j][k];
		de_3=1+2*delta_3[i][j][k];
		gam_2=1+2*gama_2[i][j][k];

		vpz=vp2;
		vpx=vp2*ep_2;
		vpy=vp2*ep_1;
		vpn1=vp2*de_1;
		vpn2=vp2*de_2;
		vpn3=vpx*de_3;//vpx*del_3
		vsz1=vs2;
		vsz2=vs2*gam_1/gam_2;
		vsz3=vs2*gam_1;
		C23_44=sqrt(vpz-vsz2)*sqrt(vpn1-vsz2);
		C12_66=sqrt(vpx-vsz3)*sqrt(vpn3-vsz3);
		C13_55=sqrt(vpz-vsz1)*sqrt(vpn2-vsz1);
		//deri calculation
		px=0;py=0;pz=0;
		qx=0;qy=0;qz=0;
		rx=0;ry=0;rz=0;
		for(l=-_m;l<=_m;l++)
		{
		    int il=i+l;
		    int jl=j+l;
		    int kl=k+l;
		    int lm=l+_m;
		    if(il>=0&&il<nypad)
		    {
			py+=coeff_2dy[lm]*p2[il][j][k];
			qy+=coeff_2dy[lm]*q2[il][j][k];
			ry+=coeff_2dy[lm]*r2[il][j][k];
		    }
		    if(jl>=0&&jl<nxpad)
		    {
			px+=coeff_2dx[lm]*p2[i][jl][k];
			qx+=coeff_2dx[lm]*q2[i][jl][k];
			rx+=coeff_2dx[lm]*r2[i][jl][k];
		    }
		    if(kl>=0&&kl<nzpad)
		    {
			pz+=coeff_2dz[lm]*p2[i][j][kl];
			qz+=coeff_2dz[lm]*q2[i][j][kl];
			rz+=coeff_2dz[lm]*r2[i][j][kl];
		    }
		}
                int im, jm, km;
                im=i+_m;
                jm=j+_m;
                km=k+_m;
		p3[im][jm][km]=2*p2[im][jm][km] - p1[im][jm][km] + dt2*(vpx*px + vsz3*py + vsz1*pz + C12_66*qx + C13_55*rx);
		q3[im][jm][km]=2*q2[im][jm][km] - q1[im][jm][km] + dt2*(C12_66*py + vsz3*qx + vpy*qy + vsz2*qz + C23_44*ry);
		r3[im][jm][km]=2*r2[im][jm][km] - r1[im][jm][km] + dt2*(C13_55*pz + C23_44*qz + vsz1*rx + vsz2*ry + vpz*rz);
	    }
}
