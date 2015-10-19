/* 3-D three-components wavefield modeling using elasic wave equation in tilted ORT media.

   Copyright (C) 2012 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng, Tengfei Wang and Wei Kang
   Modified: Yanadet Sripanich (The University of Texas at Austin)
     
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

/* prepared head files by myself */
#include "_fd.h"
#include "_cjb.h"

/* head files automatically produced from *.c */
#include "rickerjb.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"
#include "vmodel.h"

/* wavefield propagators */
#include "fwportelasticcij.h"

int   ny,nx,nz,nt;
float dx,dy,dz,dt,dxm,dym;
bool snapshot;

int main(int  argc,char **argv)
{
    int   isx,isy,isz,bd;

    int   i,j,k,im,jm,it,it1,it2,its,depth;
	int   nth, rank;
    float t,f0;
    float fx,fy,fz,dt2;

    float ***c11, ***c22, ***c33, ***c12, ***c13, ***c23, ***c44, ***c55, ***c66;
    float ***phaix, ***phaiy, ***phaiz;
    const char *source;
    bool mig;
    
    sf_init(argc,argv);
    
    // Output files
    sf_file Fo1, Fo2, Fo3;
    sf_file snapx, snapy, snapz;

    float t0=0.01;       // time delay of the wavelet(if f0=30Hz, t0=0.04s)*/
    float A=1.0;           // the amplitude of wavelet 
    if (!sf_getfloat("freq",&f0)) f0=40.0;    // main frequency of the wavelet(usually 30Hz)
    source = sf_getstring("source"); // source location
    if (NULL == (source = sf_getstring("source"))) source="m";

    clock_t t1, t2, t3;
    float   timespent;

    t1=clock();
    
    /* time samping paramter */
    if (!sf_getint("nt",&nt))   nt=301;
    if (!sf_getfloat("dt",&dt)) dt=0.001;
    if (!sf_getint("bd",&bd)) bd=20;
    if (!sf_getbool("mig",&mig)) mig=false;
    sf_warning("nt=%d dt=%f",nt,dt);
    
    /* setup I/O files */
    sf_file Fc11, Fc22, Fc33, Fc12, Fc13, Fc23, Fc44, Fc55, Fc66;
    sf_file Fphiz, Fphiy, Fphix;
    sf_file Fdatax, Fdatay, Fdataz;
    
    Fc11 = sf_input ("c11");  /* c11 using standard input */
    Fc22 = sf_input ("c22");  /* c22 */
    Fc33 = sf_input ("c33");  /* c33 */
    Fc12 = sf_input ("c12");  /* c12 */
    Fc13 = sf_input ("c13");  /* c13 */
    Fc23 = sf_input ("c23");  /* c23 */
    Fc44 = sf_input ("c44");  /* c44 */
    Fc55 = sf_input ("c55");  /* c55 */
    Fc66 = sf_input ("c66");  /* c66 */
    Fphix = sf_input ("phix");  /* phix x ccw*/
    Fphiy = sf_input ("phiy");  /* phiy y ccw*/
    Fphiz = sf_input ("phiz");  /* phiz z ccw */

    /* Read/Write axes */
    sf_axis az, ax, ay;
    az = sf_iaxa(Fc11,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fc11,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fc11,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    int nxpad, nypad, nzpad;

    nxpad=nx+2*bd;
    nypad=ny+2*bd;
    nzpad=nz+2*bd;

    sf_warning("nxpad=%d nypad=%d nzpad=%d ",nxpad,nypad,nzpad);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);
    
    float ***dataz, ***datax, ***datay;
    /* migration flag */
    if (mig) {
    Fdatax = sf_input ("datax");
    Fdatay = sf_input ("datay");
    Fdataz = sf_input ("dataz");
    if (!sf_getint("depth",&depth)) depth=0;
    int tmp;
    if (!sf_histint(Fdatax,"n1",&tmp) || nt!=tmp) sf_error("No n1=nt in datax");
    if (!sf_histint(Fdatax,"n2",&tmp) || nx!=tmp) sf_error("No n2=nx in datax");
    if (!sf_histint(Fdatax,"n3",&tmp) || ny!=tmp) sf_error("No n3=ny in datax");
    dataz = sf_floatalloc3(nt,nx,ny);
    datax = sf_floatalloc3(nt,nx,ny);
    datay = sf_floatalloc3(nt,nx,ny);
    sf_floatread(dataz[0][0],nt*nx*ny,Fdataz);
    sf_floatread(datax[0][0],nt*nx*ny,Fdatax);
    sf_floatread(datay[0][0],nt*nx*ny,Fdatay);
    } else {
    dataz = NULL;
    datax = NULL;
    datay = NULL;
    }
    
    c11=sf_floatalloc3(nzpad,nxpad,nypad);	
    c22=sf_floatalloc3(nzpad,nxpad,nypad);
    c33=sf_floatalloc3(nzpad,nxpad,nypad);	
    c12=sf_floatalloc3(nzpad,nxpad,nypad);
    c13=sf_floatalloc3(nzpad,nxpad,nypad);	
    c23=sf_floatalloc3(nzpad,nxpad,nypad);
    c44=sf_floatalloc3(nzpad,nxpad,nypad);
    c55=sf_floatalloc3(nzpad,nxpad,nypad);	
    c66=sf_floatalloc3(nzpad,nxpad,nypad);	
    phaix=sf_floatalloc3(nzpad,nxpad,nypad);
    phaiy=sf_floatalloc3(nzpad,nxpad,nypad);
    phaiz=sf_floatalloc3(nzpad,nxpad,nypad);

    /* read velocity model */
    for(i=bd;i<nypad-bd;i++)
        for(j=bd;j<nxpad-bd;j++){
          sf_floatread(&c11[i][j][bd],nz,Fc11);
          sf_floatread(&c22[i][j][bd],nz,Fc22);
          sf_floatread(&c33[i][j][bd],nz,Fc33);
          sf_floatread(&c12[i][j][bd],nz,Fc12);
          sf_floatread(&c13[i][j][bd],nz,Fc13);
          sf_floatread(&c23[i][j][bd],nz,Fc23);
          sf_floatread(&c44[i][j][bd],nz,Fc44);
          sf_floatread(&c55[i][j][bd],nz,Fc55);
          sf_floatread(&c66[i][j][bd],nz,Fc66);
          sf_floatread(&phaix[i][j][bd],nz,Fphix);
          sf_floatread(&phaiy[i][j][bd],nz,Fphiy);
          sf_floatread(&phaiz[i][j][bd],nz,Fphiz);
       }

    vmodelboundary3d(c11, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c22, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c33, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c12, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c13, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c23, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c44, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c55, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(c66, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(phaix, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(phaiy, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(phaiz, nx, ny, nz, nxpad, nypad, nzpad, bd);

    for(i=0;i<nypad;i++)
        for(j=0;j<nxpad;j++)
          for(k=0;k<nzpad;k++){
             phaix[i][j][k] *= SF_PI/180.0;
             phaiy[i][j][k] *= SF_PI/180.0;
             phaiz[i][j][k] *= SF_PI/180.0;
          }
    sf_warning("Read velocity model parameters ok !");

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

    /* Export wavefield snapshots */
    if (!sf_getbool("snapshot",&snapshot)) snapshot=false;
    int snapstep;
    if (snapshot) {
        if (!sf_getint("snapstep",&snapstep)) snapstep=1;
        
        int nsnap = (int) floor(nt/snapstep);
        
        /* setup I/O files */
        snapx = sf_output("snapx"); /* snapshots of original elasticwave iLine x-component */
        snapy = sf_output("snapy"); /* snapshots of original elasticwave iLine y-component */
        snapz = sf_output("snapz"); /* snapshots of original elasticwave xLine z-component */
        puthead4x(snapx, nz, nx, ny, nsnap, dz/1000.0, dx/1000.0, dy/1000.0, dt, 0.0, 0.0, 0.0, 0.0);
        puthead4x(snapy, nz, nx, ny, nsnap, dz/1000.0, dx/1000.0, dy/1000.0, dt, 0.0, 0.0, 0.0, 0.0);
        puthead4x(snapz, nz, nx, ny, nsnap, dz/1000.0, dx/1000.0, dy/1000.0, dt, 0.0, 0.0, 0.0, 0.0);
       
    }

    /* source definition */
     switch(source[0])  {
        case 'm':
            isy=nypad/2;
            isx=nxpad/2;
            isz=nzpad/2;
            break;
        case 'f':
            isy=nypad/2;
            isx=nxpad/2;
            isz=bd+(nzpad-2*bd)/4;
            break;
        case 's':
            isy=nypad/2;
            isx=nxpad/2;
            isz=bd;
            break;
     }

    dt2=dt*dt;

#ifdef _OPENMP
    #pragma omp parallel
	{
	  nth = omp_get_num_threads();
	  rank = omp_get_thread_num();
	  sf_warning("Using %d threads, this is %dth thread",nth, rank);
	}
#endif

    float*** px_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** pz_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** qx_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** qz_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** rx_tmp=sf_floatalloc3(nzpad,nxpad,nypad);
    float*** rz_tmp=sf_floatalloc3(nzpad,nxpad,nypad);

    if (mig) {
        it1 = nt -1;
        it2 = -1;
        its = -1;
    } else {
        it1 = 0;
        it2 = nt;
        its = 1;
    }
	/*********the kernel calculation ************/
	for(it=it1;it!=it2;it+=its)
	{
	     t=it*dt;
             
         if (mig) { // inject data
             for(i=0;i<ny;i++)
             for(j=0;j<nx;j++)
             {
                    p2[i+bd][j+bd][depth]=datax[i][j][it];
                    q2[i+bd][j+bd][depth]=datay[i][j][it];
                    r2[i+bd][j+bd][depth]=dataz[i][j][it];
             }
         } else { // inject source
         /* source Type 0: oriented 45 degree to vertical and 45 degree azimuth: Yan & Sava (2012) */
             p2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // x-component
             q2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // y-component
             r2[isy][isx][isz]+=Ricker(t, f0, t0, A);  // z-component

             // 3D exploding force source (e.g., Wu's PhD

/*               for(k=-1;k<=1;k++)*/
/*               for(i=-1;i<=1;i++)*/
/*               for(j=-1;j<=1;j++)*/
/*               {*/
/*                if(fabs(i)+fabs(j)+fabs(k)==3)*/
/*                {*/
/*                     p2[isy+k][isx+i][isz+j]+=i*Ricker(t, f0, t0, A);  // x-component*/
/*                     q2[isy+k][isx+i][isz+j]+=k*Ricker(t, f0, t0, A);  // y-component*/
/*                     r2[isy+k][isx+i][isz+j]+=j*Ricker(t, f0, t0, A);  // z-component*/
/*                }*/
/*               }*/
         }
         
         if (!mig) { // Modelling
             if (snapshot && it%snapstep == 0) { // output snapshots for all time steps
                for(i=0;i<ny;i++) {
                    im=i+bd;
        	        for(j=0;j<nx;j++) {
                        jm=j+bd;
                        sf_floatwrite(&p3[im][jm][bd],nz,snapx);
                        sf_floatwrite(&q3[im][jm][bd],nz,snapy);
                        sf_floatwrite(&r3[im][jm][bd],nz,snapz);
                        }
                }
             }
         }
         
         
         
  	     fwportelastic3d(dt2,p1,p2,p3,q1,q2,q3,r1,r2,r3,
				         px_tmp,pz_tmp,
				         qx_tmp,qz_tmp,
				         rx_tmp,rz_tmp,
                         coeff_2dx,coeff_2dy,coeff_2dz,
                         coeff_1dx,coeff_1dy,coeff_1dz,
                         dx,dy,dz,nxpad,nypad,nzpad,
			 c11,c22,c33,c12,c13,c23,c44,c55,c66,phaix,phaiy,phaiz);

         if (mig) { // Migration
             if (snapshot && it%snapstep == 0) { // output snapshots for all time steps
                for(i=0;i<ny;i++) {
                    im=i+bd;
        	        for(j=0;j<nx;j++) {
                        jm=j+bd;
                        sf_floatwrite(&p3[im][jm][bd],nz,snapx);
                        sf_floatwrite(&q3[im][jm][bd],nz,snapy);
                        sf_floatwrite(&r3[im][jm][bd],nz,snapz);
                        }
                }
             }
         }
         
         
         
            if(it==nt-1) { // output snapshot at nt-1
                // output iLine 
    	     	for(i=0;i<ny;i++) {
                    im=i+bd;
    		        for(j=0;j<nx;j++) {
                        jm=j+bd;
                        sf_floatwrite(&p3[im][jm][bd],nz,Fo1);
                        sf_floatwrite(&q3[im][jm][bd],nz,Fo2);
                        sf_floatwrite(&r3[im][jm][bd],nz,Fo3);
                        }
                }
            }
          
          // Forward to the next time step
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
           if(mig) {
               sf_warning("Backward propagation...  it= %d   t=%f",it,t);
           } else {
               sf_warning("Forward propagation...  it= %d   t=%f",it,t);
           }
     }

    printf("ok3\n");

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for 3D ORT elastic modeling: %f(second)",timespent);

    free(**p1);
    free(**p2);
    free(**p3);
    free(**q1);
    free(**q2);
    free(**q3);
    free(**r1);
    free(**r2);
    free(**r3);
    free(**px_tmp);
    free(**qx_tmp);
    free(**rx_tmp);
    free(**pz_tmp);
    free(**qz_tmp);
    free(**rz_tmp);

    free(**c11);
    free(**c33);
    free(**c13);
    free(**c55);
    free(**c66);
    free(**phaiz);
    free(**phaiy);
    free(**phaix);
		
    return 0;
}
