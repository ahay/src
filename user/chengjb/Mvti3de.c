/* 3-D three-components wavefield modeling using elasic wave equation in VTI media.

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
#include "vmodel.h"

/* wavefield propagators */
#include "fwpvtielastic.h"

int   ny,nx,nz,ns;
float dx,dy,dz,dt,dxm,dym;

int main(int  argc,char **argv)
{
    int   isx,isy,isz,bd;

    int   i,j,k,im,jm,it;
	int   nth, rank;
    float t;
    float fx,fy,fz,dt2;

    float ***vp0, ***vs0, ***epsi, ***delta, ***gama;

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
    if (!sf_getint("bd",&bd)) bd=20;

    sf_warning("ns=%d dt=%f",ns,dt);

    /* setup I/O files */
    sf_file Fvp0, Fvs0, Fep, Fde, Fga;

    Fvp0 = sf_input ("in");  /* vp0 using standard input */
    Fvs0 = sf_input ("vs0");  /* vs0 */
    Fep = sf_input ("epsi");  /* epsi */
    Fde = sf_input ("del");  /* delta */
    Fga = sf_input ("gam");  /* gama */

    /* Read/Write axes */
    sf_axis az, ax, ay;
    az = sf_iaxa(Fvp0,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fvp0,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fvp0,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    int nxpad, nypad, nzpad;

    nxpad=nx+2*bd;
    nypad=ny+2*bd;
    nzpad=nz+2*bd;

    sf_warning("nxpad=%d nypad=%d nzpad=%d ",nxpad,nypad,nzpad);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);

    vp0=sf_floatalloc3(nzpad,nxpad,nypad);	
    vs0=sf_floatalloc3(nzpad,nxpad,nypad);	
    epsi=sf_floatalloc3(nzpad,nxpad,nypad);	
    delta=sf_floatalloc3(nzpad,nxpad,nypad);	
    gama=sf_floatalloc3(nzpad,nxpad,nypad);	

    /* read velocity model */
    for(i=bd;i<nypad-bd;i++)
        for(j=bd;j<nxpad-bd;j++){
          sf_floatread(&vp0[i][j][bd],nz,Fvp0);
          sf_floatread(&vs0[i][j][bd],nz,Fvs0);
          sf_floatread(&epsi[i][j][bd],nz,Fep);
          sf_floatread(&delta[i][j][bd],nz,Fde);
          sf_floatread(&gama[i][j][bd],nz,Fga);
       }

    vmodelboundary3d(vp0, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(vs0, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(epsi, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(delta, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(gama, nx, ny, nz, nxpad, nypad, nzpad, bd);

    sf_warning("read velocity model parameters ok");

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
    float*** xtmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** ytmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** ztmp=sf_floatalloc3(nzpad,nxpad,nypad);

	for(it=0;it<ns;it++)
	{
	     t=it*dt;

		 /* source Type 0: oriented 45 degree to vertical and 45 degree azimuth: Yan & Sava (2012) */
         p2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // x-component
         q2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // y-component
         r2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // z-component

         // 3D exploding force source (e.g., Wu's PhD)
		 // Displacement source equivalent to point-source in pseudo-pure-mode wavefield
		 /*
         for(k=-1;k<=1;k++)
            for(i=-1;i<=1;i++)
               for(j=-1;j<=1;j++)
               {
				// source Type 1
                if(fabs(i)+fabs(j)+fabs(k)==3)
                     p2[isy+k][isx+i][isz+j]+=i*Ricker(t, f0, t0, A);  // x-component
                     q2[isy+k][isx+i][isz+j]+=k*Ricker(t, f0, t0, A);  // y-component
                     r2[isy+k][isx+i][isz+j]+=j*Ricker(t, f0, t0, A);  // z-component
                }
				// source Type 2
                if(i+j+k==3||i+j+k==-3)
                {
                     p2[isy+k][isx+i][isz+j]+=i*Ricker(t, f0, t0, A);  // x-component
                     q2[isy+k][isx+i][isz+j]+=k*Ricker(t, f0, t0, A);  // y-component
                     r2[isy+k][isx+i][isz+j]+=j*Ricker(t, f0, t0, A);  // z-component
                }
               }
		 */
  	     fwpvtielastic3d(dt2,p1,p2,p3,q1,q2,q3,r1,r2,r3,xtmp,ytmp,ztmp,
                         coeff_2dx,coeff_2dy,coeff_2dz,
                         coeff_1dx,coeff_1dy,coeff_1dz,
                         dx,dy,dz,nxpad,nypad,nzpad,
	                     vp0,vs0,epsi,delta,gama);

         if(it==ns-1) // output snapshot
         {
	     	for(i=0;i<ny;i++)
                {
                    im=i+bd;
		            for(j=0;j<nx;j++)
                    {
                        jm=j+bd;
                        sf_floatwrite(&p3[im][jm][bd],nz,Fo1);
                        sf_floatwrite(&q3[im][jm][bd],nz,Fo2);
                        sf_floatwrite(&r3[im][jm][bd],nz,Fo3);
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

    free(**vp0);
    free(**vs0);
    free(**epsi);
    free(**gama);
    free(**delta);
		
    return 0;
}
