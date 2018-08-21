/* 3-D three-components wavefield modeling using pure mode SH-wave equation in 3D VTI media.
 * Note: The z-components of pure-mode qSV-wave are zero.

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
#include "vmodel.h"

/* wave-mode separation operators */
/* wavefield propagators */
#include "fwpvtipseudosh.h"

static int   ny,nx,nz,ns;
static float dx,dy,dz,dt;

int main(int  argc,char **argv)
{
    int   isx,isy,isz;

    int   i,j,k,im,jm,it,ns,bd;
    float t;
    float dx,dz,dy,dt;

    float ***vp0, ***vs0, ***epsi, ***del, ***gam;

    sf_init(argc,argv);

    sf_file Fo1, Fo2, Fo3;

    float f0=40;         // main frequency of the wavelet(usually 30Hz)
    float t0=0.04;       // time delay of the wavelet(if f0=30Hz, t0=0.04s)*/
    float A=1.0;           // the amplitude of wavelet

    clock_t t2, t3;
    float   timespent;

    /* time samping paramter */
    if (!sf_getint("ns",&ns)) ns=301;
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

    sf_warning("_m=%d _mix=%d",_m,_mix);
    sf_warning("nx=%d ny=%d nz=%d ",nx,ny,nz);
    sf_warning("dx=%f dy=%f dz=%f ",dx,dy,dz);

    sf_warning("bd=%d",bd);

    int nxpad, nypad, nzpad;

    nxpad=nx+2*bd;
    nypad=ny+2*bd;
    nzpad=nz+2*bd;

    vp0=sf_floatalloc3(nzpad,nxpad,nypad);
    vs0=sf_floatalloc3(nzpad,nxpad,nypad);
    epsi=sf_floatalloc3(nzpad,nxpad,nypad);
    del=sf_floatalloc3(nzpad,nxpad,nypad);
    gam=sf_floatalloc3(nzpad,nxpad,nypad);

    /* read velocity model */
    for(i=bd;i<nypad-bd;i++)
        for(j=bd;j<nxpad-bd;j++){
          sf_floatread(&vp0[i][j][bd],nz,Fvp0);
          sf_floatread(&vs0[i][j][bd],nz,Fvs0);
          sf_floatread(&epsi[i][j][bd],nz,Fep);
          sf_floatread(&del[i][j][bd],nz,Fde);
          sf_floatread(&gam[i][j][bd],nz,Fga);
       }

    vmodelboundary3d(vp0, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(vs0, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(epsi, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(del, nx, ny, nz, nxpad, nypad, nzpad, bd);
    vmodelboundary3d(gam, nx, ny, nz, nxpad, nypad, nzpad, bd);

    sf_warning("read velocity model parameters ok");

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
	    float*** p1=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** p2=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** p3=sf_floatalloc3(nzpad,nxpad,nypad);

        zero3float(p1,nzpad,nxpad,nypad);
        zero3float(p2,nzpad,nxpad,nypad);
        zero3float(p3,nzpad,nxpad,nypad);

        float*** r1=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** r2=sf_floatalloc3(nzpad,nxpad,nypad);
        float*** r3=sf_floatalloc3(nzpad,nxpad,nypad);

        zero3float(r1,nzpad,nxpad,nypad);
        zero3float(r2,nzpad,nxpad,nypad);
        zero3float(r3,nzpad,nxpad,nypad);

        t2=clock();

        /* setup I/O files */
        Fo1 = sf_output("out"); /* pure SH-wave x-component */
        Fo2 = sf_output("SHy"); /* pure SH-wave z-component */
        Fo3 = sf_output("SH"); /* scalar SH-wave */

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
	     r2[isy][isx][isz]+=Ricker(t, f0, t0, A);
		 sf_warning("p2=%f r2=%f\n", p2[isy][isx][isz],r2[isy][isx][isz]);

  	     fwpvti3dpseudosh(dt,p1,p2,p3,r1,r2,r3,
                          coeff_2dx,coeff_2dy,coeff_2dz,coeff_1dx,coeff_1dy,coeff_1dz,
	                      vp0,vs0,epsi,del,gam,
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
                        sf_floatwrite(&r3[im][jm][bd],nz,Fo2);
						for(k=bd;k<nz+bd;k++)
						{
						   float tmp=p3[im][jm][k]+r3[im][jm][k];
						   sf_floatwrite(&tmp, 1, Fo3);
						}
                    }
                }
             }
            for(i=0;i<nypad;i++)
            for(j=0;j<nxpad;j++)
            for(k=0;k<nzpad;k++)
            {
                    p1[i][j][k]=p2[i][j][k];
                    p2[i][j][k]=p3[i][j][k];

                    r1[i][j][k]=r2[i][j][k];
                    r2[i][j][k]=r3[i][j][k];
           }
           sf_warning("forward propagation...  it= %d",it);
     }

    printf("ok3\n");

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for 3D VTI pseudo-pure modeling: %f(second)",timespent);

    free(**p1);
    free(**p2);
    free(**p3);
    free(**r1);
    free(**r2);
    free(**r3);

    free(**vp0);
    free(**vs0);
    free(**epsi);
    free(**del);
    free(**gam);

    return 0;
}
