/* 3-D three-components wavefield modeling using pseudo-pure mode P-wave equation in tilted ORT media.

   Refernces:
             Cheng et al. (15th IWSA, 2012);
             Cheng and Kang (SEG Abstract, 2012);
             Kang and Cheng (SEG Abstract, 2012)
             Wang et al.(SEG Abstract, 2012)      

   Copyright (C) 2012 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng, Tengfei Wang and Wei Kang
     
   This code is first written by Tengfei Wang at Tongji University,
   and then optimzied by Jiubing Cheng for Madagascar version at BEG,
   University of Texas at Austin.

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

/* head files aumatically produced from *.c */
#include "ricker.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"

/* wave-mode separation operators */
/* wavefield propagators */
#include "fwportpseudop.h"

int   ny,nx,nz,ns;
float dx,dy,dz,dt,dxm,dym;

int main(int  argc,char **argv)
{
    int   isx,isy,isz;

    int   i,j,k,im,jm,it,ns,bd;
    float t;
    float fx,fy,fz,dx,dz,dy,dt;

    float vp0, vs0, epsi1, epsi2, del1, del2, del3, gam1, gam2;

    sf_init(argc,argv);

    sf_file Fo1, Fo2, Fo3;

    float f0=40;         // main frequency of the wavelet(usually 30Hz)
    float t0=0.04;       // time delay of the wavelet(if f0=30Hz, t0=0.04s)*/
    float A=1.0;           // the amplitude of wavelet 

    clock_t t1, t2, t3;
    float   timespent;

    t1=clock();

    /* time samping paramter */
    if (!sf_getint("ns",&ns)) ns=301;
    if (!sf_getint("ny",&ny)) ny=101;
    if (!sf_getint("nx",&nx)) nx=101;
    if (!sf_getint("nz",&nz)) nz=101;
    if (!sf_getfloat("dt",&dt)) dt=0.001;
    if (!sf_getfloat("dx",&dx)) dx=0.0;
    if (!sf_getfloat("dy",&dy)) dy=0.0;
    if (!sf_getfloat("dz",&dz)) dz=0.0;

    if (!sf_getfloat("vp0",&vp0)) vp0=3000.0;
    if (!sf_getfloat("vs0",&vs0)) vs0=1500.0;
    if (!sf_getfloat("epsi1",&epsi1)) epsi1=0.2;
    if (!sf_getfloat("epsi2",&epsi2)) epsi2=0.067;
    if (!sf_getfloat("del1",&del1)) del1=0.1;
    if (!sf_getfloat("del2",&del2)) del2=-0.0422;
    if (!sf_getfloat("del3",&del3)) del3=0.125;
    if (!sf_getfloat("gam1",&gam1)) gam1=0.1;
    if (!sf_getfloat("gam2",&gam2)) gam2=0.047;

    if (!sf_getint("bd",&bd)) bd=20;

    sf_warning("ns=%d dt=%f",ns,dt);

    fy=0.0;
    fx=0.0;
    fz=0.0;

    sf_warning("_m=%d _mix=%d",_m,_mix);
    sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);
    sf_warning("vp0=%f vs0=%f ",vp0, vs0);
    sf_warning("epsi1=%f epsi2=%f ",epsi1, epsi2);
    sf_warning("gam1=%f gam2=%f ",gam1, gam2);
    sf_warning("del1=%f del2=%f del3=%f",del1, del2, del3);

    sf_warning("bd=%d",bd);

    int mm=2*_m+1;
    int mmix=2*_mix+1;

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

        int nxpad, nypad, nzpad;

        nxpad=nx+2*bd;
        nypad=ny+2*bd;
        nzpad=nz+2*bd;

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
        Fo1 = sf_output("out");          /* pseudo-pure P-wave iLine x-component */
        Fo2 = sf_output("PseudoPurePy"); /* pseudo-pure P-wave iLine y-component */
        Fo3 = sf_output("PseudoPurePz"); /* pseudo-pure P-wave xLine z-component */

        puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
        puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
        puthead3x(Fo3, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

 /*****************************************************************************
 *  Calculating polarization deviation operator for wave-mode separation
 * ***************************************************************************/
        /* source definition */
        isy=nypad/2;
        isx=nxpad/2;
        isz=nzpad/2;

	/*********the kernel calculation ************/
	for(it=0;it<ns;it++)
	{
	     t=it*dt;
             //sf_warning("t=%f f0=%f t0=%f A=%f",t, f0, t0, A);
	     p2[isy][isx][isz]+=Ricker(t, f0, t0, A);
	     q2[isy][isx][isz]+=Ricker(t, f0, t0, A);
	     r2[isy][isx][isz]+=Ricker(t, f0, t0, A);
		 sf_warning("p2=%f q2=%f r2=%f\n", p2[isy][isx][isz],q2[isy][isx][isz],r2[isy][isx][isz]);

  	     fwportpseudophomo(dt,p1,p2,p3,q1,q2,q3,r1,r2,r3,
                           coeff_2dx,coeff_2dy,coeff_2dz,coeff_1dx,coeff_1dy,coeff_1dz,
	                       vp0,vs0,epsi1,del1,gam1,epsi2,del2,gam2,del3,
                           nxpad,nypad,nzpad,dx,dy,dz);//for ORT

        if(it==ns-1) // output snapshot
        {
        // output iLine 
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

    return 0;
}
