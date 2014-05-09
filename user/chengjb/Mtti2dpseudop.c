/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in TTI media.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang
     
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

/* prepared head files by myself */
#include "_fd.h"
#include "_cjb.h"

/* head files aumatically produced from *.c */
#include "ricker.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"
#include "kykxkz2yxz.h"
#include "kykxkztaper.h"
#include "clipsmthspec.h"

/* wave-mode separation operators */
#include "devttip.h"
#include "filter2dsep.h"

/*wavefield propagators */
#include "fwpttipseudop.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
    int	ix, iz, jx, jz, ixf, izf, ixx, izz, i,j,im, jm,nx,nz,nxf,nzf,nxpad,nzpad,it,ii,jj;
    float   kxmax,kzmax;

    float   A, f0, t, t0, dx, dz, dxf, dzf, dt, dkx, dkz, dt2, div;
    int     mm, nvx, nvz, ns;
    int     hnkx, hnkz, nkx, nkz, nxz, nkxz;
    int     hnkx1=1, hnkz1=1, nkx1, nkz1;
    int     isx, isz, isxm, iszm; /*source location */
    int     itaper; /* tapering or not for spectrum of oprtator*/

    int     nstep;            /* every nstep in spatial grids to calculate filters sparsely*/

    float   *coeff_1dx, *coeff_1dz, *coeff_2dx, *coeff_2dz; /* finite-difference coefficient */

    float **apvx, **apvz, **apvxx, **apvzz;    /* projection deviation operator of P-wave for a location */

    float ****ex=NULL, ****ez=NULL;                      /* operator for whole model for P-wave*/
    float **exx=NULL, **ezz=NULL;                        /* operator for constant model for P-wave*/

    float **vp0, **vs0, **epsi, **del, **theta;         /* velocity model */
    float **p1, **p2, **p3, **q1, **q2, **q3, **p3c=NULL, **q3c=NULL, **sum=NULL;  /* wavefield array */
    float *kx, *kz, *kkx, *kkz, *kx2, *kz2, **taper;

    clock_t t2, t3, t4, t5;
    float   timespent, fx,fz; 
    char    *tapertype;

    int     isep=1;
    int     ihomo=1;

    double  vp2, vs2, ep2, de2, the;

    sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6, Fo7, Fo8;
    sf_file Fvp0, Fvs0, Feps, Fdel, Fthe;

    sf_axis az, ax;
       
    sf_init(argc,argv);

    /* t1=clock(); */
 
    /*  wavelet parameter for source definition */
    f0=30.0;                  
    t0=0.04;                  
    A=1.0;                  

    /* time samping paramter */
    if (!sf_getint("ns",&ns)) ns=301;
    if (!sf_getfloat("dt",&dt)) dt=0.001;
    if (!sf_getint("isep",&isep)) isep=0;             /* if isep=1, separate wave-modes */
    if (!sf_getint("ihomo",&ihomo)) ihomo=0;          /* if ihomo=1, homogeneous medium */
    if (NULL== (tapertype=sf_getstring("tapertype"))) tapertype="D"; /* taper type*/
    if (!sf_getint("nstep",&nstep)) nstep=1; /* grid step to calculate operators: 1<=nstep<=5 */

    sf_warning("isep=%d",isep);
    sf_warning("ihomo=%d",ihomo);
    sf_warning("tapertype=%s",tapertype);
    sf_warning("nstep=%d",nstep);

    sf_warning("ns=%d dt=%f",ns,dt);
    sf_warning("read velocity model parameters");

    /* setup I/O files */
    Fvp0 = sf_input ("in");  /* vp0 using standard input */
    Fvs0 = sf_input ("vs0");  /* vs0 */
    Feps = sf_input ("epsi");  /* epsi */
    Fdel = sf_input ("del");  /* delta */
    Fthe = sf_input ("the");  /* theta */

    /* Read/Write axes */
    az = sf_iaxa(Fvp0,1); nvz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fvp0,2); nvx = sf_n(ax); dx = sf_d(ax)*1000.0;
    fx=sf_o(ax);
    fz=sf_o(az);

    /* source definition */
    isx=nvx/2;
    isz=nvz/2;
    /* isz=nvz*2/5; */

    /* wave modeling space */
    nx=nvx;
    nz=nvz;
    nxpad=nx+2*_m;
    nzpad=nz+2*_m;

    sf_warning("dx=%f dz=%f",dx,dz);

    sf_warning("nx=%d nz=%d nxpad=%d nzpad=%d", nx,nz,nxpad,nzpad);

    vp0=sf_floatalloc2(nz,nx);	
    vs0=sf_floatalloc2(nz,nx);	
    epsi=sf_floatalloc2(nz,nx);	
    del=sf_floatalloc2(nz,nx);	
    theta=sf_floatalloc2(nz,nx);	

    nxz=nx*nz;
    mm=2*_m+1;

    dt2=dt*dt;
    isxm=isx+_m;  /* source's x location */
    iszm=isz+_m;  /* source's z-location */

    /* read velocity model */
    sf_floatread(vp0[0],nxz,Fvp0);
    sf_floatread(vs0[0],nxz,Fvs0);
    sf_floatread(epsi[0],nxz,Feps);
    sf_floatread(del[0],nxz,Fdel);
    sf_floatread(theta[0],nxz,Fthe);

    for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
	    theta[i][j] *= SF_PI/180.0;

    t2=clock();

    /* setup I/O files */
    Fo1 = sf_output("out"); /* pseudo-pure P-wave x-component */
    Fo2 = sf_output("PseudoPurePz"); /* pseudo-pure P-wave z-component */
    Fo3 = sf_output("PseudoPureP"); /* scalar P-wave field using divergence operator */
    puthead3(Fo1, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz, fx, 0.0);
    puthead3(Fo2, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz, fx, 0.0);
    puthead3(Fo3, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz, fx, 0.0);

    /*****************************************************************************
     *  Calculating polarization deviation operator for wave-mode separation
     * ***************************************************************************/
    if(isep==1)
    {
	/* calculate spatial steps for operater in sparsely sampling grid point */
	dxf=dx*nstep;
	dzf=dz*nstep;
	nxf=nx/nstep+1;
	nzf=nz/nstep+1;

	/* operators length for calculation */
	hnkx=400.0/dx;
	hnkz=400.0/dz;
	nkx=2*hnkx+1;   /* operator length in kx-direction */
	nkz=2*hnkz+1;   /* operator length in kz-direction */

	/* truncated spatial operators length for filtering*/
	hnkx1=200.0/dx;
	hnkz1=200.0/dz;
	nkx1=2*hnkx1+1;
	nkz1=2*hnkz1+1;

	sf_warning("nx=%d nz=%d nxf=%d nzf=%d", nx,nz,nxf,nzf);
	sf_warning("dx=%f dz=%f dxf=%f dzf=%f", dx,dz,dxf,dzf);

	sf_warning("hnkx=%d hnkz=%d nkx=%d nkz=%d", hnkx, hnkz, nkx, nkz);
	sf_warning("hnkx1=%d hnkz1=%d nkx1=%d nkz1=%d", hnkx1, hnkz1, nkx1, nkz1);

	dkx=2*SF_PI/dx/nkx;
	dkz=2*SF_PI/dz/nkz;
	kxmax=SF_PI/dx;
	kzmax=SF_PI/dz;

	kx=sf_floatalloc(nkx);
	kz=sf_floatalloc(nkx);
	kkx=sf_floatalloc(nkx);
	kkz=sf_floatalloc(nkx);
	kx2=sf_floatalloc(nkx);
	kz2=sf_floatalloc(nkx);

	taper=sf_floatalloc2(nkz, nkx);

	/* define axis samples and taper in wavenumber domain */
	kxkztaper(kx, kz, kkx, kkz, kx2, kz2, taper, nkx, nkz, hnkx, hnkz, dkx, dkz, kxmax, kzmax, tapertype);
           
	/* truncation of spatial filter */
	p3c=sf_floatalloc2(nz,nx);
	q3c=sf_floatalloc2(nz,nx);
	sum=sf_floatalloc2(nz,nx);

	if(ihomo==1)
	{
	    exx=sf_floatalloc2(nkz1, nkx1);
	    ezz=sf_floatalloc2(nkz1, nkx1);
	}
	else{ /* to store spatail varaied operators */
	    ex=sf_floatalloc4(nkz1, nkx1, nzf, nxf);
	    ez=sf_floatalloc4(nkz1, nkx1, nzf, nxf);
	}
	nkxz=nkx*nkz;

	apvx=sf_floatalloc2(nkz, nkx);
	apvz=sf_floatalloc2(nkz, nkx);
	apvxx=sf_floatalloc2(nkz, nkx);
	apvzz=sf_floatalloc2(nkz, nkx);

	/* setup I/O files */
	Fo4 = sf_output("apvx"); /* P-wave projection deviation x-comp */
	Fo5 = sf_output("apvz"); /* P-wave projection deviation z-comp */
	Fo6 = sf_output("apvxx"); /* P-wave projection deviation x-comp in (x,z) domain */
	Fo7 = sf_output("apvzz"); /* P-wave projection deviation z-comp in (x,z) domain */

	puthead1(Fo4, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
	puthead1(Fo5, nkz, nkx, dkz, -kzmax, dkx, -kxmax);

	puthead2(Fo6, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);
	puthead2(Fo7, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);

	/*************calculate projection deviation grid-point by grid-point **********/
	for(ix=0,ixf=0;ix<nx;ix+=nstep,ixf++)
	{
	    for(iz=0,izf=0;iz<nz;iz+=nstep,izf++)
	    {
	        vp2=vp0[ix][iz]*vp0[ix][iz];
	        vs2=vs0[ix][iz]*vs0[ix][iz];
	        ep2=1.0+2*epsi[ix][iz];
	        de2=1.0+2*del[ix][iz];
	        the=theta[ix][iz];

                if(ixf%10==0&&izf%100==0) sf_warning("Operator: nxf=%d ixf=%d izf=%d vp2=%f vs2=%f",nxf, ixf,izf,vp2,vs2);

   	        /*************calculate projection deviation without tapering **********/
                itaper=0;
                /* devvtip: projection deviation operators for P-wave in TTI media */
                devttip(apvx,apvz,kx,kz,kkx,kkz,taper,hnkx,hnkz,dkx,dkz,vp2,vs2,ep2,de2,the,itaper);

                /* inverse Fourier transform */
                kxkz2xz(apvx, apvxx, hnkx, hnkz, nkx, nkz);
                kxkz2xz(apvz, apvzz, hnkx, hnkz, nkx, nkz);

                /* truncation and saving of operator in space-domain */
                if(ihomo==1)
                {
		    for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++)
			for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++)
			{
			    exx[jx+hnkx1][jz+hnkz1]=apvxx[ixx][izz];
			    ezz[jx+hnkx1][jz+hnkz1]=apvzz[ixx][izz];
			}
                }else{
		    for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++)
			for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++)
			{
			    ex[ixf][izf][jx+hnkx1][jz+hnkz1]=apvxx[ixx][izz];
			    ez[ixf][izf][jx+hnkx1][jz+hnkz1]=apvzz[ixx][izz];
			}
                }
                if((ix==nx/2&&iz==nz/2&&ihomo==0)||ihomo==1)
                {
		    sf_floatwrite(apvx[0], nkxz, Fo4);
		    sf_floatwrite(apvz[0], nkxz, Fo5);

		    sf_floatwrite(apvxx[0], nkxz, Fo6);
		    sf_floatwrite(apvzz[0], nkxz, Fo7);
                }
                if(ihomo==1) goto loop;

	    }/* iz loop */
	}/* ix loop */
    loop:;

	free(kx);
	free(kz);
	free(kx2);
	free(kz2);
	free(kkx);
	free(kkz);

	free(*taper);

	free(*apvx);
	free(*apvz);
	free(*apvxx);
	free(*apvzz);
    }/* isep */
    /****************End of Calculating Projection Deviation Operator****************/
    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("Computation time (operators): %f (second)",timespent);

    /****************begin to calculate wavefield****************/
    coeff_2dx=sf_floatalloc(mm);
    coeff_2dz=sf_floatalloc(mm);
    coeff_1dx=sf_floatalloc(mm);
    coeff_1dz=sf_floatalloc(mm);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dz,dz);
    coeff1d(coeff_1dx,dx);
    coeff1d(coeff_1dz,dz);

    p1=sf_floatalloc2(nzpad, nxpad);
    p2=sf_floatalloc2(nzpad, nxpad);
    p3=sf_floatalloc2(nzpad, nxpad);

    q1=sf_floatalloc2(nzpad, nxpad);
    q2=sf_floatalloc2(nzpad, nxpad);
    q3=sf_floatalloc2(nzpad, nxpad);

    zero2float(p1, nzpad, nxpad);
    zero2float(p2, nzpad, nxpad);
    zero2float(p3, nzpad, nxpad);
        
    zero2float(q1, nzpad, nxpad);
    zero2float(q2, nzpad, nxpad);
    zero2float(q3, nzpad, nxpad);
        
    if(isep==1)
    {
        /* setup I/O files */
        Fo8 = sf_output("PseudoPureSepP"); /* scalar P-wave field using polarization projection oprtator*/

        puthead3(Fo8, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz, fx, 0.0);
    } else {
	Fo8 = NULL;
    }

    sf_warning("==================================================");
    sf_warning("==  Propagation Using Pseudo-Pure P-Wave Eq.    ==");
    sf_warning("==================================================");
    t4=clock();
    for(it=0;it<ns;it++)
    {
	t=it*dt;

	p2[isxm][iszm]+=Ricker(t, f0, t0, A);
	q2[isxm][iszm]+=Ricker(t, f0, t0, A);

	/* fwpttipseudop: forward-propagating in TTI media with pseudo-pure P-wave equation */
	fwpttipseudop(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz,
		      dx, dz, nx, nz, nxpad, nzpad, vp0, vs0, epsi, del, theta);

	/******* output wavefields: component and divergence *******/
	if(it==ns-1)
	{
	    for(i=0;i<nx;i++)
	    {
		im=i+_m;
		for(j=0;j<nz;j++)
		{
		    jm=j+_m;
		    sf_floatwrite(&p3[im][jm],1,Fo1);
		    sf_floatwrite(&q3[im][jm],1,Fo2);

		    div=p3[im][jm]+q3[im][jm];

		    sf_floatwrite(&div,1,Fo3);
		}
	    }/* i loop*/

	    /* correct projection error for accurate separate qP wave in spatial domain */
	    if(isep==1)
	    {
		zero2float(p3c,nz,nx);
		zero2float(q3c,nz,nx);
		zero2float(sum, nz, nx);

		if(ihomo==1)
		    filter2dsepglobal(p3, q3, p3c, q3c, exx, ezz, nx, nz, hnkx1, hnkz1);
		else
		    filter2dsep(p3, q3, p3c, q3c, ex, ez, nx, nz, nstep, hnkx1, hnkz1);

		for(i=0;i<nx;i++)
		    for(j=0;j<nz;j++)
			sum[i][j]=p3c[i][j]+q3c[i][j];

		sf_floatwrite(sum[0],nx*nz, Fo8);
	    }
	}/* (it+1)%ntstep==0 */

	/**************************************/
	for(i=0,ii=_m;i<nx;i++,ii++)
	    for(j=0,jj=_m;j<nz;j++,jj++)
	    {
		p1[ii][jj]=p2[ii][jj];	
		p2[ii][jj]=p3[ii][jj];	

		q1[ii][jj]=q2[ii][jj];	
		q2[ii][jj]=q3[ii][jj];	
	    }

	if(it%50==0)
	    sf_warning("Pseudo: it= %d",it);
    }/* it loop */
    t5=clock();
    timespent=(float)(t5-t4)/CLOCKS_PER_SEC;
    sf_warning("Computation time (propagation + separation): %f (second)",timespent);

    if(isep==1)
    {
	free(*p3c);
	free(*q3c);
	free(*sum);

	if(ihomo==1)
	{
	    free(*exx);
	    free(*ezz);
	}else{
	    free(***ex);
	    free(***ez);
	}
    }

    free(*p1);
    free(*p2);
    free(*p3);
    free(*q1);
    free(*q2);
    free(*q3);

    free(*vp0);
    free(*vs0);
    free(*epsi);
    free(*del);
    free(*theta);

    exit(0);
}
