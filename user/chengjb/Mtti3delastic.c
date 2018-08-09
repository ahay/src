/* 3-D three-components wavefield modeling using elasic wave equation in tilted TI media.

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
#include <time.h>
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* prepared head files by myself */
#include "_cjb.h"
#include "_fd.h"

/* head files aumatically produced from *.c */
#include "ricker.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"

/* wavefield propagators */
#include "fwpttielastic.h"

static int   ny,nx,nz,ns;
static float dx,dy,dz,dt;

int main(int  argc,char **argv)
{
    int   isx,isy,isz;

    int   i,j,k,im,jm,it;
    float t;
    float fx,fy,fz,dt2;

    float ***vp0, ***vs0, ***epsi, ***delta, ***gama;
    float ***theta, ***phai;

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

    sf_warning("ns=%d dt=%f",ns,dt);

    /* setup I/O files */
    sf_file Fvp0, Fvs0, Fep, Fde, Fga;
    sf_file Fthe, Fphi;

    Fvp0 = sf_input ("in");  /* vp0 using standard input */
    Fvs0 = sf_input ("vs0");  /* vs0 */
    Fep = sf_input ("epsi");  /* epsi */
    Fde = sf_input ("delta");  /* delta */
    Fga = sf_input ("gama");  /* gama */
    Fthe = sf_input ("theta");  /* theta */
    Fphi = sf_input ("phai");  /* phi */

    /* Read/Write axes */
    sf_axis az, ax, ay;
    az = sf_iaxa(Fvp0,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fvp0,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fvp0,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);

    vp0=sf_floatalloc3(nz,nx,ny);	
    vs0=sf_floatalloc3(nz,nx,ny);	
    epsi=sf_floatalloc3(nz,nx,ny);	
    delta=sf_floatalloc3(nz,nx,ny);	
    gama=sf_floatalloc3(nz,nx,ny);	
    theta=sf_floatalloc3(nz,nx,ny);
    phai=sf_floatalloc3(nz,nx,ny);

    /* read velocity model */
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          sf_floatread(vp0[i][j],nz,Fvp0);
          sf_floatread(vs0[i][j],nz,Fvs0);
          sf_floatread(epsi[i][j],nz,Fep);
          sf_floatread(delta[i][j],nz,Fde);
          sf_floatread(gama[i][j],nz,Fga);
          sf_floatread(theta[i][j],nz,Fthe);
          sf_floatread(phai[i][j],nz,Fphi);

          for(k=0;k<nz;k++){
             theta[i][j][k] *= PI/180.0;
             phai[i][j][k] *= PI/180.0;
          }
        }
    sf_warning("read velocity model parameters ok");

    int mm=2*m+1;
    int mmix=2*mix+1;
 
    sf_warning("m=%d mix=%d",m,mix);

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

        //calc_coeff_1d(coeff_1dx);
        //calc_coeff_1d(coeff_1dy);
        //calc_coeff_1d(coeff_1dz);

        int nxpad, nypad, nzpad;

        nxpad=nx+2*m;
        nypad=ny+2*m;
        nzpad=nz+2*m;

	float*** p1=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** p2=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** p3=sf_floatalloc3(nzpad,nxpad,nypad);

        zero3float(p1,nzpad,nxpad,nypad);
        zero3float(p2,nzpad,nxpad,nypad);
        zero3float(p3,nzpad,nxpad,nypad);
    
        float*** q1=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** q2=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** q3=sf_floatalloc3(nzpad,nxpad,nypad);

        zero3float(q1,nzpad,nxpad,nypad);
        zero3float(q2,nzpad,nxpad,nypad);
        zero3float(q3,nzpad,nxpad,nypad);

        float*** r1=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** r2=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** r3=sf_floatalloc3(nzpad,nxpad,nypad);

        zero3float(r1,nzpad,nxpad,nypad);
        zero3float(r2,nzpad,nxpad,nypad);
        zero3float(r3,nzpad,nxpad,nypad);

        t2=clock();

        /* setup I/O files */
        Fo1 = sf_output("out");      /* original elasticwave iLine x-component */
        Fo2 = sf_output("Elasticy"); /* original elasticwave iLine y-component */
        Fo3 = sf_output("Elasticz"); /* original elasticwave xLine z-component */

        puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
        puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
        puthead3x(Fo3, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

        /* source definition */
        isy=nypad/2;
        isx=nxpad/2;
        isz=nzpad/2;

        dt2=dt*dt;

	/*********the kernel calculation ************/
	for(it=0;it<ns;it++)
	{
	     t=it*dt;
             //sf_warning("t=%f f0=%f t0=%f A=%f",t, f0, t0, A);
             
             // 3D exploding force source (e.g., Wu's PhD
             if(fabs(t-t0)<0.002) {
               for(k=-1;k<=1;k++)
               for(i=-1;i<=1;i++)
               for(j=-1;j<=1;j++)
               {
                if(fabs(i)+fabs(j)+fabs(k)==3)
                {
                     //p2[isy+k][isx+i][isz+j]+=i*Ricker(t, f0, t0, A);  // x-component
                     //q2[isy+k][isx+i][isz+j]+=k*Ricker(t, f0, t0, A);  // y-component
                     //r2[isy+k][isx+i][isz+j]+=j*Ricker(t, f0, t0, A);  // z-component
                     p2[isy+k][isx+i][isz+j]+=i;  // x-component
                     q2[isy+k][isx+i][isz+j]+=k;  // y-component
                     r2[isy+k][isx+i][isz+j]+=j;  // z-component
                }
                sf_warning("p2=%f",p2[isy+k][isx+i][isz+j]);
                sf_warning("q2=%f",q2[isy+k][isx+i][isz+j]);
                sf_warning("r2=%f",r2[isy+k][isx+i][isz+j]);
               }
             }
            
             sf_warning("=================================================t=%f ",t);
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
  	     fwpttielastic3d(dt2,p1,p2,p3,q1,q2,q3,r1,r2,r3,
                             coeff_2dx,coeff_2dy,coeff_2dz,
                             coeff_1dx,coeff_1dy,coeff_1dz,
                             dx,dy,dz,nx,ny,nz,nxpad,nypad,nzpad,
	                     vp0,vs0,epsi,delta,gama,theta,phai);

             if(it==ns-1) // output snapshot
             {
                // output iLine 
		for(i=0;i<ny;i++)
		//for(i=ny/2;i<=ny/2;i++)
                {
                    im=i+m;
		    for(j=0;j<nx;j++)
                    {
                        jm=j+m;
                        sf_floatwrite(&p3[im][jm][m],nz,Fo1);
                        sf_floatwrite(&q3[im][jm][m],nz,Fo2);
                        sf_floatwrite(&r3[im][jm][m],nz,Fo3);
                    }
                }
             }
            for(i=0;i<nypad;i++)
            for(j=0;j<nxpad;j++)
            for(k=0;k<nzpad;k++)
            {
                    p1[i][j][k]=p2[i][j][k];
                    p2[i][j][k]=p3[i][j][k];

                    q1[i][j][k]=q2[i][j][k];
                    q2[i][j][k]=q3[i][j][k];

                    r1[i][j][k]=r2[i][j][k];
                    r2[i][j][k]=r3[i][j][k];
           }

           sf_warning("forward propagation...  it= %d",it);
     }

    printf("ok3\n");

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for 3D ORT pseudo-pure modeling: %f(second)",timespent);

    free(**p1);
    free(**p2);
    free(**p3);
    free(**q1);
    free(**q2);
    free(**q3);
    free(**r1);
    free(**r2);
    free(**r3);

    free(**vp0);
    free(**vs0);
    free(**epsi);
    free(**gama);
    free(**delta);
    free(**theta);
    free(**phai);
		
    return 0;
}
