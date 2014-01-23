/* 3-D three-components wavefield modeling using elasic wave equation in TTI media.

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
#include <omp.h>

/* prepared head files by myself */
#include "_fd.h"
#include "_cjb.h"

/* head files aumatically produced from *.c */
#include "ricker.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"

/* wavefield propagators */
#include "fwpttielastic.h"

int   ny,nx,nz,ns;
float dx,dy,dz,dt;

int main(int  argc,char **argv)
{
    int   isx,isy,isz;

    int   i,j,k,it;
	int   nth, rank;
    float t;
    float fx,fy,fz,dt2;

    float vp0, vs0, epsi, delta, gama, theta, phai;

    sf_init(argc,argv);

    sf_file Fo1, Fo2, Fo3;

    float f0=40;         // main frequency of the wavelet(usually 30Hz)
    float t0=0.04;       // time delay of the wavelet(if f0=30Hz, t0=0.04s)*/
    float A=1.0;           // the amplitude of wavelet 

    clock_t t1, t2, t3;
    float   timespent;

    t1=clock();

    /* time samping paramter */
    if (!sf_getint("ns",&ns))   ns=301;
    if (!sf_getfloat("dt",&dt)) dt=0.001;
    if (!sf_getfloat("vp0",&vp0)) vp0=3000.0;
    if (!sf_getfloat("vs0",&vs0)) vs0=2000.0;
    if (!sf_getfloat("epsi",&epsi)) epsi=0.1;
    if (!sf_getfloat("delta",&delta)) delta=0.05;
    if (!sf_getfloat("gama",&gama)) gama=0.05;
    if (!sf_getfloat("theta",&theta)) theta=30.0;
    if (!sf_getfloat("phai",&phai)) phai=45.0;
    if (!sf_getint("nx",&nx))   nx=201;
    if (!sf_getint("ny",&ny))   ny=201;
    if (!sf_getint("nz",&nz))   nz=201;
    if (!sf_getfloat("dx",&dx)) dx=5.0;
    if (!sf_getfloat("dy",&dy)) dy=5.0;
    if (!sf_getfloat("dz",&dz)) dz=5.0;

    fy=0.0;
    fx=0.0;
    fz=0.0;

	theta *= SF_PI/180.0;
	phai *= SF_PI/180.0;

    sf_warning("ns=%d dt=%f",ns,dt);
    sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);
	sf_warning("vp0=%f vs0=%f epsi=%f delta=%f gama=%f theta=%f phai=%f",vp0, vs0, epsi, delta, gama, theta, phai);

    int mm=2*_m+1;
    int mmix=2*_mix+1;
 
    sf_warning("m=%d mix=%d",_m,_mix);

    float *coeff_2dx,*coeff_2dy,*coeff_2dz,*coeff_1dx,*coeff_1dy,*coeff_1dz;

    coeff_2dy=sf_floatalloc(mm);
    coeff_2dx=sf_floatalloc(mm);
    coeff_2dz=sf_floatalloc(mm);
    coeff_1dy=sf_floatalloc(mmix);
    coeff_1dx=sf_floatalloc(mmix);
    coeff_1dz=sf_floatalloc(mmix);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dy,dy);
    coeff2d(coeff_2dz,dz);

    coeff1dmix(coeff_1dx, dx);
    coeff1dmix(coeff_1dy, dy);
    coeff1dmix(coeff_1dz, dz);

	float*** p1=sf_floatalloc3(nz,nx,ny);
    float*** p2=sf_floatalloc3(nz,nx,ny);
    float*** p3=sf_floatalloc3(nz,nx,ny);

    zero3float(p1,nz,nx,ny);
    zero3float(p2,nz,nx,ny);
    zero3float(p3,nz,nx,ny);
    
    float*** q1=sf_floatalloc3(nz,nx,ny);
    float*** q2=sf_floatalloc3(nz,nx,ny);
    float*** q3=sf_floatalloc3(nz,nx,ny);

    zero3float(q1,nz,nx,ny);
    zero3float(q2,nz,nx,ny);
    zero3float(q3,nz,nx,ny);

    float*** r1=sf_floatalloc3(nz,nx,ny);
    float*** r2=sf_floatalloc3(nz,nx,ny);
    float*** r3=sf_floatalloc3(nz,nx,ny);

    zero3float(r1,nz,nx,ny);
    zero3float(r2,nz,nx,ny);
    zero3float(r3,nz,nx,ny);

    t2=clock();

    /* setup I/O files */
    Fo1 = sf_output("out");      /* original elasticwave iLine x-component */
    Fo2 = sf_output("Elasticy"); /* original elasticwave iLine y-component */
    Fo3 = sf_output("Elasticz"); /* original elasticwave xLine z-component */

    puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    puthead3x(Fo3, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

    /* source definition */
    isy=ny/2;
    isx=nx/2;
    isz=nz/2;

    dt2=dt*dt;

    #pragma omp parallel
	{
	  nth = omp_get_num_threads();
	  rank = omp_get_thread_num();
	  sf_warning("Using %d threads, this is %dth thread",nth, rank);
	}
    float***xtmp=sf_floatalloc3(nz,nx,ny);
    float***ytmp=sf_floatalloc3(nz,nx,ny);
    float***ztmp=sf_floatalloc3(nz,nx,ny);

	/*********the kernel calculation ************/
	for(it=0;it<ns;it++)
	{
	     t=it*dt;
             //sf_warning("t=%f f0=%f t0=%f A=%f",t, f0, t0, A);
             
             // 3D exploding force source (e.g., Wu's PhD
               for(k=-1;k<=1;k++)
               for(i=-1;i<=1;i++)
               for(j=-1;j<=1;j++)
               {
                if(fabs(i)+fabs(j)+fabs(k)==3)
                {
                     p2[isy+k][isx+i][isz+j]+=i*Ricker(t, f0, t0, A);  // x-component
                     q2[isy+k][isx+i][isz+j]+=k*Ricker(t, f0, t0, A);  // y-component
                     r2[isy+k][isx+i][isz+j]+=j*Ricker(t, f0, t0, A);  // z-component
                }
               }
             // 3D equil-energy force source (e.g., Wu's PhD)
             /*
             for(k=-1;k<=1;k++)
             for(i=-1;i<=1;i++)
             for(j=-1;j<=1;j++)
             {
                if(fabs(i)+fabs(j)+fabs(k)==3)
                {
                   // need complement 
                }
             }
             */
  	     fwpttielastic3dhomo(dt2,p1,p2,p3,q1,q2,q3,r1,r2,r3,xtmp,ytmp,ztmp,
                         coeff_2dx,coeff_2dy,coeff_2dz,
                         coeff_1dx,coeff_1dy,coeff_1dz,
                         dx,dy,dz,nx,ny,nz,
	                     vp0,vs0,epsi,delta,gama,theta,phai);

         if(it==ns-1) // output snapshot
         {
            // output iLine 
	     	for(i=0;i<ny;i++)
                {
		            for(j=0;j<nx;j++)
                    {
                        sf_floatwrite(&p3[i][j][0],nz,Fo1);
                        sf_floatwrite(&q3[i][j][0],nz,Fo2);
                        sf_floatwrite(&r3[i][j][0],nz,Fo3);
                    }
                }
             }
            for(i=0;i<ny;i++)
            for(j=0;j<nx;j++)
            for(k=0;k<nz;k++)
            {
                    p1[i][j][k]=p2[i][j][k];
                    p2[i][j][k]=p3[i][j][k];

                    q1[i][j][k]=q2[i][j][k];
                    q2[i][j][k]=q3[i][j][k];

                    r1[i][j][k]=r2[i][j][k];
                    r2[i][j][k]=r3[i][j][k];
           }

           sf_warning("forward propagation...  it= %d   t=%f",it,t);
     }

    printf("ok3\n");

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for 3D TTI elastic modeling: %f(second)",timespent);

    free(**p1);
    free(**p2);
    free(**p3);
    free(**q1);
    free(**q2);
    free(**q3);
    free(**r1);
    free(**r2);
    free(**r3);
    free(**xtmp);
    free(**ytmp);
    free(**ztmp);

    return 0;
}
