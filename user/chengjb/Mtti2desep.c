/* 2-D two-components wavefield modeling using original elastic displacement wave equation in TTI media.

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
#include "kykxkztaper.h"
#include "kykxkz2yxz.h"
#include "clipsmthspec.h"

/* wave-mode separation operators */
#include "polttipsv.h"
#include "filter2dsep.h"

/* wavefield propagators */
#include "fwpttielastic.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
    int	ix, iz, jx, jz,ixx,izz,ixf,izf,i,j,im, jm,nx,nz,nxf,nzf,nxpad,nzpad,it,ii,jj;
    float   kxmax,kzmax;

    float   f0, t, t0, dx, dz, dxf, dzf,dt, dkx, dkz, dt2;
    int     mm, nvx, nvz, ns;
    int     hnkx, hnkz, nkx, nkz, nxz, nkxz;
    int     hnkx1=0, hnkz1=0, nkx1, nkz1;
    int     isx, isz, isxm, iszm; /*source location */

    int     itaper; /* tapering or not for spectrum of oprtator*/
    int     nstep;            /* every nstep in spatial grids to calculate filters sparsely*/

    float   *coeff_1dx, *coeff_1dz, *coeff_2dx, *coeff_2dz; /* finite-difference coefficient */

    float **apx, **apz, **apxx, **apzz;        /* polarization operator of P-wave for a location */
    float **apxs, **apzs, **apxxs, **apzzs;    /* polarization operator of SV-wave for a location */

    float ****ex=NULL, ****ez=NULL;                      /* operator for whole model for P-wave*/
    float ****exs=NULL, ****ezs=NULL;                    /* operator for whole model for SV-wave*/
    float **exx=NULL, **ezz=NULL;                        /* operator for constant model for P-wave*/
    float **exxs=NULL, **ezzs=NULL;                      /* operator for constant model for SV-wave*/

    float **vp0, **vs0, **epsi, **del, **theta;         /* velocity model */
    float **p1, **p2, **p3, **q1, **q2, **q3, **p3c, **q3c, **sum;  /* wavefield array */

    float *kx, *kz, *kkx, *kkz, *kx2, *kz2, **taper;

    float   A, fx, fz;

    int     isep=1;
    int     ihomo=1;

    char    *tapertype;

    double  vp2, vs2, ep2, de2, the;

    sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6, Fo7, Fo8, Fo9, Fo10, Fo11, Fo12;
    sf_file Fvp0, Fvs0, Feps, Fdel, Fthe;

    sf_axis az, ax;

    sf_init(argc,argv);
   
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
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    /* source definition */
    isx=nvx/2;
    isz=nvz/2;
    /* isz=nvz*2/5; */

    /* wave modeling space */
    nx=nvx;
    nz=nvz;
    nxpad=nx+2*_m;
    nzpad=nz+2*_m;

    sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);

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

    Fo1 = sf_output("out"); /* Elastic-wave x-component */
    Fo2 = sf_output("Elasticz"); /* Elastic-wave z-component */
    /* setup I/O files */
    puthead3(Fo1, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0);
    puthead3(Fo2, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0);

    /*****************************************************************************
     *  Calculating polarization operator for wave-mode separation
     * ***************************************************************************/
    if(isep==1)
    {
	sf_warning("==================================================");
	sf_warning("==      Calculating Polarization Operator       ==");
	sf_warning("==================================================");
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

	nkxz=nkx*nkz;

	/* truncation of spatial filter */
	if(ihomo==1)
	{
	    exx=sf_floatalloc2(nkz1, nkx1);
	    ezz=sf_floatalloc2(nkz1, nkx1);
	    exxs=sf_floatalloc2(nkz1, nkx1);
	    ezzs=sf_floatalloc2(nkz1, nkx1);
	}else{
	    ex=sf_floatalloc4(nkz1, nkx1, nz, nx);
	    ez=sf_floatalloc4(nkz1, nkx1, nz, nx);
	    exs=sf_floatalloc4(nkz1, nkx1, nz, nx);
	    ezs=sf_floatalloc4(nkz1, nkx1, nz, nx);
	}
	/*****************************************************************************
	 *  Calculating polarization operator for wave-mode separation
	 * ***************************************************************************/
	apx=sf_floatalloc2(nkz, nkx);
	apz=sf_floatalloc2(nkz, nkx);

	apxs=sf_floatalloc2(nkz, nkx);
	apzs=sf_floatalloc2(nkz, nkx);

	apxx=sf_floatalloc2(nkz, nkx);
	apzz=sf_floatalloc2(nkz, nkx);

	apxxs=sf_floatalloc2(nkz, nkx);
	apzzs=sf_floatalloc2(nkz, nkx);

	/* setup I/O files for wavenumber-domain operators */
	Fo3 = sf_output("apx"); /*  P-wave's polarization x-comp */
	Fo4 = sf_output("apz"); /*  P-wave's polarization z-comp */
	Fo5 = sf_output("apxs"); /*  SV-wave's polarization x-comp */
	Fo6 = sf_output("apzs"); /*  SV-wave's polarization z-comp */

	puthead1(Fo3, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
	puthead1(Fo4, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
	puthead1(Fo5, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
	puthead1(Fo6, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
 
	/* setup I/O files for space-domain operators */
	Fo7 = sf_output("apxx");  /* P-wave's polarization x-comp in (x,z) domain */
	Fo8 = sf_output("apzz");  /* P-wave's polarization z-comp in (x,z) domain */
	Fo9 = sf_output("apxxs"); /* SV-wave's polarization x-comp in (x,z) domain */
	Fo10 = sf_output("apzzs"); /* SV-wave's polarization z-comp in (x,z) domain */

	puthead2(Fo7, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);
	puthead2(Fo8, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);
	puthead2(Fo9, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);
	puthead2(Fo10, nkz, nkx, dz/1000.0, 0.0, dx/1000.0, 0.0);

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

   	        /*************calculate projection operrate with tapering **********/
                zero2float(apx, nkz, nkx);
                zero2float(apz, nkz, nkx);
                zero2float(apxs, nkz, nkx);
                zero2float(apzs, nkz, nkx);
               
                /* polvtipsv: P- and SV-wave polarization operators in VTI media  */
                itaper=1;
                polttipsv(apx,apz,apxs,apzs,kx,kz,kkx,kkz,taper,hnkx,hnkz,dkx,dkz,
                          vp2,vs2,ep2,de2,the,itaper);

                ikxkz2xz(apx, apxx, hnkx, hnkz, nkx, nkz);
                ikxkz2xz(apz, apzz, hnkx, hnkz, nkx, nkz);
                ikxkz2xz(apxs, apxxs, hnkx, hnkz, nkx, nkz);
                ikxkz2xz(apzs, apzzs, hnkx, hnkz, nkx, nkz);

                /* truncation and saving of operator in space-domain */
                if(ihomo==1)
                {
		    for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++) 
			for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++) 
			{
			    exx[jx+hnkx1][jz+hnkz1]=apxx[ixx][izz]; 
			    ezz[jx+hnkx1][jz+hnkz1]=apzz[ixx][izz]; 
			    exxs[jx+hnkx1][jz+hnkz1]=apxxs[ixx][izz]; 
			    ezzs[jx+hnkx1][jz+hnkz1]=apzzs[ixx][izz]; 
			}
                }else{
		    for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++) 
			for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++) 
			{
			    ex[ixf][izf][jx+hnkx1][jz+hnkz1]=apxx[ixx][izz]; 
			    ez[ixf][izf][jx+hnkx1][jz+hnkz1]=apzz[ixx][izz]; 
			    exs[ixf][izf][jx+hnkx1][jz+hnkz1]=apxxs[ixx][izz]; 
			    ezs[ixf][izf][jx+hnkx1][jz+hnkz1]=apzzs[ixx][izz]; 
			}
                }
                
                if((ixf==nxf/2&&izf==nzf/2&&ihomo==0)||ihomo==1)
                {
		    /* write-disk operators in kx-kz domain */
		    sf_floatwrite(apx[0], nkxz, Fo3);
		    sf_floatwrite(apz[0], nkxz, Fo4);
		    sf_floatwrite(apxs[0], nkxz, Fo5);
		    sf_floatwrite(apzs[0], nkxz, Fo6);

		    /* write-disk operators in x-z domain */
		    sf_floatwrite(apxx[0], nkxz, Fo7);
		    sf_floatwrite(apzz[0], nkxz, Fo8);
		    sf_floatwrite(apxxs[0], nkxz, Fo9);
		    sf_floatwrite(apzzs[0], nkxz, Fo10);
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

	free(*apx);
	free(*apz);
	free(*apxs);
	free(*apzs);
	free(*apxx);
	free(*apzz);
	free(*apxxs);
	free(*apzzs);
    }/* isep loop */
    /****************End of Calculating Projection Deviation Operator****************/
 
    /****************begin to calculate wavefield****************/
    /****************begin to calculate wavefield****************/
    sf_warning("==================================================");
    sf_warning("==      Propagation Using Elastic Wave Eq.      ==");
    sf_warning("==================================================");
    coeff_2dx=sf_floatalloc(mm);
    coeff_2dz=sf_floatalloc(mm);
    coeff_1dx=sf_floatalloc(mm);
    coeff_1dz=sf_floatalloc(mm);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dz,dz);

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
        
    coeff1dmix(coeff_1dx,dx);
    coeff1dmix(coeff_1dz,dz);

    if(isep==1)
    {
	Fo11 = sf_output("ElasticSepP"); /*  scalar wavefield using P-wave's polarization projection oprtator*/
	Fo12 = sf_output("ElasticSepSV"); /*  scalar wavefield using SV-wave's polarization projection oprtator*/

	puthead3(Fo11, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0);
	puthead3(Fo12, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, 0.0);

	p3c=sf_floatalloc2(nz,nx);
	q3c=sf_floatalloc2(nz,nx);
	sum=sf_floatalloc2(nz,nx);

    } else {

	Fo11 = NULL;
	Fo12 = NULL;
	
	p3c=NULL;
	q3c=NULL;
	sum=NULL;

    }

    for(it=0;it<ns;it++)
    {
	t=it*dt;

	/* 2D exploding force source (e.g., Wu's PhD) */
	for(i=-1;i<=1;i++)
	    for(j=-1;j<=1;j++)
	    {
		if(SF_ABS(i)+SF_ABS(j)==2)
		{
		    p2[isxm+i][iszm+j]+=i*Ricker(t, f0, t0, A);
		    q2[isxm+i][iszm+j]+=j*Ricker(t, f0, t0, A);
		}
	    }
        /* 2D equil-energy force source (e.g., Wu's PhD) */
        /*
	  for(i=-1;i<=1;i++)
	  for(j=-1;j<=1;j++)
	  {
	  if(fabs(i)+fabs(j)==2)
	  {
	  if(i==-1&&j==1)  
	  q2[isxm+i][iszm+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	  if(i==-1&&j==-1) 
	  p2[isxm+i][iszm+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	  if(i==1&&j==1)  
	  p2[isxm+i][iszm+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	  if(i==1&&j==-1) 
	  q2[isxm+i][iszm+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	  }
	  }
        */
	/* fwpvtielastic: forward-propagating using original elastic equation of displacement in VTI media*/
	fwpttielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
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
		}
	    }/* i loop*/

                       
	    if(isep==1)
	    {
		/* applying P-wave polarization projection operator in spatial domain */
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

		sf_floatwrite(sum[0],nx*nz, Fo11);
          
		/* applying SV-wave polarization projection operator in spatial domain */
		zero2float(p3c,nz,nx);
		zero2float(q3c,nz,nx);
		zero2float(sum, nz, nx);

		if(ihomo==1)
		    filter2dsepglobal(p3, q3, p3c, q3c, exxs, ezzs, nx, nz, hnkx1, hnkz1);
		else
		    filter2dsep(p3, q3, p3c, q3c, exs, ezs, nx, nz, nstep, hnkx1, hnkz1);

		for(i=0;i<nx;i++)
		    for(j=0;j<nz;j++)
			sum[i][j]=p3c[i][j]+q3c[i][j];

		sf_floatwrite(sum[0],nx*nz, Fo12);
	    }/* isep==1 */

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

	if(it%100==0)
	    sf_warning("Elastic: it= %d",it);
    }/* it loop */

    if(isep==1)
    {
	free(*p3c);
	free(*q3c);
	free(*sum);

	if(ihomo==1)
	{
	    free(*exx);
	    free(*ezz);
	    free(*exxs);
	    free(*ezzs);
	}else{
	    free(***ex);
	    free(***ez);
	    free(***exs);
	    free(***ezs);
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
