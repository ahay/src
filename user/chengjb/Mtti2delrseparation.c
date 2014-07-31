/* 2-D two-components elastic displacement wavefield snapshot modeling  
   and wave mode separation based on low-rank approximate mixed-domain
   separators for TTI media.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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
#include <assert.h>

/* prepared head files by myself */
#include "_cjb.h"

#ifndef M
#define M 5   /* 10th-order finite-difference: accuracy is 2*M */
#endif
#ifndef Mix
#define Mix 5   /* order of finite-difference for mix-derivative (2*mix) */
#endif

/* head files aumatically produced from C programs */
#include "zero.h"
#include "ricker.h"
#include "puthead.h"
#include "kykxkztaper.h"
#include "fdcoef.h"
#include "fwpttielastic.h"
#include "seplowrank.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
   sf_init(argc,argv);

   clock_t t1, t2, t3, t4, t5;
   float   timespent;

   int   ns, it;
   float dt;

   t1=clock();

   if (!sf_getint("ns",&ns)) ns=301;   // total time sample
   if (!sf_getfloat("dt",&dt)) dt=0.001;   // t-axis sampling step during wavefield extrapolation

   sf_warning("ns=%d dt=%f",ns,dt);

   sf_warning("##### Step 1: Read anisotropic models``");
   /* setup I files */
   sf_file Fvp, Fvs, Fep, Fde, Fth;
   sf_axis az, ax;

   Fvp = sf_input("in");  /* vp0 using standard input */
   Fvs = sf_input("vs0");
   Fep = sf_input("epsi");
   Fde = sf_input("del");
   Fth = sf_input("the");

   /* Read/Write axes */
   int nxv, nzv;
   float fx, fz;
   float dx, dz;
   az = sf_iaxa(Fvp,1); nzv = sf_n(az); dz = sf_d(az)*1000.0;
   ax = sf_iaxa(Fvp,2); nxv = sf_n(ax); dx = sf_d(ax)*1000.0;
   fx=sf_o(ax)*1000.0;
   fz=sf_o(az)*1000.0;

   /* wave modeling space */
   int nx, nz, nxz;
   nx=nxv;
   nz=nzv;
   nxz=nx*nz;

   float **vp, **vs, **ep, **de, **th;
   vp = sf_floatalloc2(nzv,nxv);
   vs = sf_floatalloc2(nzv,nxv);
   ep = sf_floatalloc2(nzv,nxv);
   de = sf_floatalloc2(nzv,nxv);
   th = sf_floatalloc2(nzv,nxv);

   /* read velocity model */
   sf_floatread(vp[0],nxz,Fvp);
   sf_floatread(vs[0],nxz,Fvs);
   sf_floatread(ep[0],nxz,Fep);
   sf_floatread(de[0],nxz,Fde);
   sf_floatread(th[0],nxz,Fth);

   int i, j;
   for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
         th[i][j] *= SF_PI/180.0;

   sf_warning("##### Step 2: Read low-rank separation matrix");
   /* Fourier spectra demension */
   int nkz,nkx,nk;
   nkx=nx;
   nkz=nz;
   /*  read low-rank decomppsed matrix for wave mode separation */
   sf_file SepPBx, SepPAx, SepPCx, SepPBz, SepPAz, SepPCz;
   sf_file SepSBx, SepSAx, SepSCx, SepSBz, SepSAz, SepSCz;
   sf_axis a1, a2;

   /////////////////////////////////////////////
   /////// low-rank decomposed matrix for P-wave
   int   m2xp, n2xp;
   float *lmatxp, *mmatxp, *rmatxp;

   ///////////////////////////////// SepPBx
   SepPBx = sf_input("SepPBx");
   a1 = sf_iaxa(SepPBx, 1);
   a2 = sf_iaxa(SepPBx, 2);
   m2xp = sf_n(a1);
   if(nxz != sf_n(a2)){
       sf_warning("SepPBx: n2 != nxz");
       exit(0);
   }
    sf_warning("SepPBx n1=m2xp=%d n2=nxz=%d ",m2xp,nxz);
    lmatxp = sf_floatalloc(m2xp*nxz);
    sf_floatread(lmatxp, m2xp*nxz, SepPBx);
    
    ///////////////////////////////// SepPAx
    SepPAx = sf_input("SepPAx");
    a1 = sf_iaxa(SepPAx, 1);
    a2 = sf_iaxa(SepPAx, 2);
    n2xp = sf_n(a1);
    sf_warning("SepPAx n1=n2xp=%d n2=%d ",n2xp,sf_n(a2));
    if(sf_n(a2) != m2xp){ 
        sf_warning("SepPAx: n2 != m2xp ");
        exit(0);
    }
    mmatxp = sf_floatalloc(n2xp*m2xp);
    sf_floatread(mmatxp, n2xp*m2xp, SepPAx);

    ///////////////////////////////// SepPCx
    SepPCx = sf_input("SepPCx");
    a1 = sf_iaxa(SepPCx, 1);
    a2 = sf_iaxa(SepPCx, 2);
    nk = sf_n(a1);
    if(nk != nkx*nkz){
        sf_warning("SepPCx: n1 != nkx*nkz");
        exit(0);
    }
    sf_warning("SepPCx n1=nk=%d n2=%d ",nk,sf_n(a2));
    if(sf_n(a2) != n2xp){ 
        sf_warning("SepPCx: n2 != n2xp ");
        exit(0);
    }
    rmatxp = sf_floatalloc(nk*n2xp);
    sf_floatread(rmatxp, nk*n2xp, SepPCx);

    int   m2zp, n2zp;
    float *lmatzp, *mmatzp, *rmatzp;

    ///////////////////////////////// SepPBz
    SepPBz = sf_input("SepPBz");
    a1 = sf_iaxa(SepPBz, 1);
    a2 = sf_iaxa(SepPBz, 2);
    m2zp = sf_n(a1);
    if(nxz != sf_n(a2)){
        sf_warning("SepPBz: n2 != nxz");
        exit(0);
    }
    sf_warning("SepPBz n1=m2zp=%d n2=nxz=%d ",m2zp,nxz);
    lmatzp = sf_floatalloc(m2zp*nxz);
    sf_floatread(lmatzp, m2zp*nxz, SepPBz);
    
    ///////////////////////////////// SepPAz
    SepPAz = sf_input("SepPAz");
    a1 = sf_iaxa(SepPAz, 1);
    a2 = sf_iaxa(SepPAz, 2);
    n2zp = sf_n(a1);
    sf_warning("SepPAz n1=n2zp=%d n2=%d ",n2zp,sf_n(a2));
    if(sf_n(a2) != m2zp){ 
        sf_warning("SepPAz: n2 != m2zp ");
        exit(0);
    }
    mmatzp = sf_floatalloc(n2zp*m2zp);
    sf_floatread(mmatzp, n2zp*m2zp, SepPAz);

    ///////////////////////////////// SepPCz
    SepPCz = sf_input("SepPCz");
    a1 = sf_iaxa(SepPCz, 1);
    a2 = sf_iaxa(SepPCz, 2);
    nk = sf_n(a1);
    if(nk != nkx*nkz){
        sf_warning("SepPCz: n1 != nkx*nkz");
        exit(0);
    }
    sf_warning("SepPCz n1=nk=%d n2=%d ",nk,sf_n(a2));
    if(sf_n(a2) != n2zp){ 
        sf_warning("SepPCz: n2 != n2zp ");
        exit(0);
    }
    rmatzp = sf_floatalloc(nk*n2zp);
    sf_floatread(rmatzp, nk*n2zp, SepPCz);

    /////////////////////////////////////////////
    /////// low-rank decomposed matrix for SV-wave
    int   m2xs, n2xs;
    float *lmatxs, *mmatxs, *rmatxs;

    ///////////////////////////////// SepSBx
    SepSBx = sf_input("SepSBx");
    a1 = sf_iaxa(SepSBx, 1);
    a2 = sf_iaxa(SepSBx, 2);
    m2xs = sf_n(a1);
    if(nxz != sf_n(a2)){
        sf_warning("SepSBx: n2 != nxz");
        exit(0);
    }
    sf_warning("SepSBx n1=m2xs=%d n2=nxz=%d ",m2xs,nxz);
    lmatxs = sf_floatalloc(m2xs*nxz);
    sf_floatread(lmatxs, m2xs*nxz, SepSBx);
    
    ///////////////////////////////// SepSAx
    SepSAx = sf_input("SepSAx");
    a1 = sf_iaxa(SepSAx, 1);
    a2 = sf_iaxa(SepSAx, 2);
    n2xs = sf_n(a1);
    sf_warning("SepSAx n1=n2xs=%d n2=%d ",n2xs,sf_n(a2));
    if(sf_n(a2) != m2xs){ 
        sf_warning("SepSAx: n2 != m2xs ");
        exit(0);
    }
    mmatxs = sf_floatalloc(n2xs*m2xs);
    sf_floatread(mmatxs, n2xs*m2xs, SepSAx);

    ///////////////////////////////// SepSCx
    SepSCx = sf_input("SepSCx");
    a1 = sf_iaxa(SepSCx, 1);
    a2 = sf_iaxa(SepSCx, 2);
    nk = sf_n(a1);
    if(nk != nkx*nkz){
        sf_warning("SepSCx: n1 != nkx*nkz");
        exit(0);
    }
    sf_warning("SepSCx n1=nk=%d n2=%d ",nk,sf_n(a2));
    if(sf_n(a2) != n2xs){ 
        sf_warning("SepSCx: n2 != n2xs ");
        exit(0);
    }
    rmatxs = sf_floatalloc(nk*n2xs);
    sf_floatread(rmatxs, nk*n2xs, SepSCx);

    int   m2zs, n2zs;
    float *lmatzs, *mmatzs, *rmatzs;

    ///////////////////////////////// SepSBz
    SepSBz = sf_input("SepSBz");
    a1 = sf_iaxa(SepSBz, 1);
    a2 = sf_iaxa(SepSBz, 2);
    m2zs = sf_n(a1);
    if(nxz != sf_n(a2)){
        sf_warning("SepSBz: n2 != nxz");
        exit(0);
    }
    sf_warning("SepSBz n1=m2zs=%d n2=nxz=%d ",m2zs,nxz);
    lmatzs = sf_floatalloc(m2zs*nxz);
    sf_floatread(lmatzs, m2zs*nxz, SepSBz);
    
    ///////////////////////////////// SepSAz
    SepSAz = sf_input("SepSAz");
    a1 = sf_iaxa(SepSAz, 1);
    a2 = sf_iaxa(SepSAz, 2);
    n2zs = sf_n(a1);
    sf_warning("SepSAz n1=n2zs=%d n2=%d ",n2zs,sf_n(a2));
    if(sf_n(a2) != m2zs){ 
        sf_warning("SepSAz: n2 != m2zs ");
        exit(0);
    }
    mmatzs = sf_floatalloc(n2zs*m2zs);
    sf_floatread(mmatzs, n2zs*m2zs, SepSAz);

    ///////////////////////////////// SepSCz
    SepSCz = sf_input("SepSCz");
    a1 = sf_iaxa(SepSCz, 1);
    a2 = sf_iaxa(SepSCz, 2);
    nk = sf_n(a1);
    if(nk != nkx*nkz){
        sf_warning("SepSCz: n1 != nkx*nkz");
        exit(0);
    }
    sf_warning("SepSCz n1=nk=%d n2=%d ",nk,sf_n(a2));
    if(sf_n(a2) != n2zs){ 
        sf_warning("SepSCz: n2 != n2zs ");
        exit(0);
    }
    rmatzs = sf_floatalloc(nk*n2zs);
    sf_floatread(rmatzs, nk*n2zs, SepSCz);

    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("CPU time for read separators: %f(second)",timespent);

   sf_warning("##### Step 3: Wavefield extrap & sep.");
   /****************begin to calculate wavefield****************/
   /*  wavelet parameter for source definition */
   float A, f0, t0, ft=0.0, dt2;
   dt2=dt*dt;
   f0=30.0;                  
   t0=0.04;                  
   A=1;                  

   int nxpad, nzpad;
   int mm=2*M+1;

   nxpad=nx+2*M;
   nzpad=nz+2*M;

   sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);
   sf_warning("nx=%d nz=%d nxpad=%d nzpad=%d", nx,nz,nxpad,nzpad);

   /* source definition */
   int ixs, izs, ixms, izms;
   ixs=nxv/2;
   izs=nzv/2;
   ixms=ixs+M;  /* source's x location */
   izms=izs+M;  /* source's z-location */

   /* setup I/O files */
   sf_file Fex, Fez, Fp, Fs;
   Fex = sf_output("out");
   Fez = sf_output("Elasticz");
   Fp  = sf_output("ElasticSepP");
   Fs  = sf_output("ElasticSepSV");

   puthead3(Fex, nz, nx, 1, dz/1000, dx/1000, dt, fz/1000, fx/1000, ft);
   puthead3(Fez, nz, nx, 1, dz/1000, dx/1000, dt, fz/1000, fx/1000, ft);
   puthead3(Fp, nz, nx, 1, dz/1000, dx/1000, dt, fz/1000, fx/1000, ft);
   puthead3(Fs, nz, nx, 1, dz/1000, dx/1000, dt, fz/1000, fx/1000, ft);

    float *coeff_2dx=sf_floatalloc(mm);
    float *coeff_2dz=sf_floatalloc(mm);
    float *coeff_1dx=sf_floatalloc(mm);
    float *coeff_1dz=sf_floatalloc(mm);

    coeff2d(coeff_2dx,dx);
    coeff2d(coeff_2dz,dz);
    coeff1dmix(coeff_1dx,dx);
    coeff1dmix(coeff_1dz,dz);

    float **p1=sf_floatalloc2(nzpad, nxpad);
    float **p2=sf_floatalloc2(nzpad, nxpad);
    float **p3=sf_floatalloc2(nzpad, nxpad);

    float **q1=sf_floatalloc2(nzpad, nxpad);
    float **q2=sf_floatalloc2(nzpad, nxpad);
    float **q3=sf_floatalloc2(nzpad, nxpad);

    zero2float(p1, nzpad, nxpad);
    zero2float(p2, nzpad, nxpad);
    zero2float(p3, nzpad, nxpad);

    zero2float(q1, nzpad, nxpad);
    zero2float(q2, nzpad, nxpad);
    zero2float(q3, nzpad, nxpad);

    sf_warning("==================================================");
    sf_warning("==  Propagation Using Elastic anisotropic Eq.   ==");
    sf_warning("==================================================");

    float *pp, *qq;
    pp=sf_floatalloc(nxz);
    qq=sf_floatalloc(nxz);

    int k, ii, jj, im, jm;
    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

    int iflag=0;

    for(it=0;it<ns;it++)
    {
	float t=it*dt;

	if(it%50==0)
		sf_warning("Elastic: it= %d",it);

		/*
	    p2[ixms][izms]+=Ricker(t, f0, t0, A);
		q2[ixms][izms]+=Ricker(t, f0, t0, A);
		*/

        // 2D exploding force source (e.g., Wu's PhD
        for(i=-1;i<=1;i++)
        for(j=-1;j<=1;j++)
        {
             if(fabs(i)+fabs(j)==2)
             {
                  p2[ixms+i][izms+j]+=i*Ricker(t, f0, t0, A);
                  q2[ixms+i][izms+j]+=j*Ricker(t, f0, t0, A);
             }
        }
        // 2D equil-energy force source (e.g., Wu's PhD)
        /*
        for(i=-1;i<=1;i++)
        for(j=-1;j<=1;j++)
        {
             if(fabs(i)+fabs(j)==2)
             {
                  if(i==-1&&j==1)  
                    q2[ixms+i][izms+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==-1&&j==-1) 
                   p2[ixms+i][izms+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==1&&j==1)  
                   p2[ixms+i][izms+j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
                  if(i==1&&j==-1) 
                    q2[ixms+i][izms+j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
             }
        }
        */

        /* fwpttielastic: forward-propagating using original elastic equation of displacement in VTI media*/
        fwpttielastic(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz, coeff_1dx, coeff_1dz,
                      dx, dz, nx, nz, nxpad, nzpad, vp, vs, ep, de, th);

        /******* output wavefields: component and divergence *******/
        if(it==ns-1)
	{
              t3=clock();

              k=0;
	      for(i=0;i<nx;i++)
              {
                   im=i+M;
		   for(j=0;j<nz;j++)
		   {
                       jm=j+M;

                       pp[k] = p3[im][jm];
                       qq[k] = q3[im][jm];

                       k++;      
		    }
              }// i loop
              sf_floatwrite(pp, nxz, Fex);
              sf_floatwrite(qq, nxz, Fez);

              /* separate qP wave  */
              sf_warning("========================================================"); 
              sf_warning("separate qP-wave based on lowrank decomp."); 
              seplowrank2d(lmatxp,rmatxp,mmatxp,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xp,n2xp,iflag);
              seplowrank2d(lmatzp,rmatzp,mmatzp,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zp,n2zp,iflag);

	      for(k=0;k<nxz;k++) pp[k] += qq[k];
              sf_floatwrite(pp, nxz, Fp);
          
              k=0;
	      for(i=0;i<nx;i++)
              {
                   im=i+M;
		   for(j=0;j<nz;j++)
		   {
                       jm=j+M;

                       pp[k] = p3[im][jm];
                       qq[k] = q3[im][jm];

                       k++;      
		    }
              }// i loop

              /* separate qSV wave  */
              sf_warning("========================================================"); 
              sf_warning("separate qSV-wave based on lowrank decomp."); 
              seplowrank2d(lmatxs,rmatxs,mmatxs,pp,ijkx,ijkz,nx,nz,nxz,nk,m2xs,n2xs,iflag);
              seplowrank2d(lmatzs,rmatzs,mmatzs,qq,ijkx,ijkz,nx,nz,nxz,nk,m2zs,n2zs,iflag);

	      for(k=0;k<nxz;k++) pp[k] += qq[k];
              sf_floatwrite(pp, nxz, Fs);

              t4=clock();
              timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
              sf_warning("CPU time for separation per timestep: %f(second)",timespent);
         }/* (it+1)%ntstep==0 */

         /**************************************/
 	 for(i=0,ii=M;i<nx;i++,ii++)
	    for(j=0,jj=M;j<nz;j++,jj++)
	    {
		p1[ii][jj]=p2[ii][jj];	
		p2[ii][jj]=p3[ii][jj];	

		q1[ii][jj]=q2[ii][jj];	
		q2[ii][jj]=q3[ii][jj];	
	    }

    }/* it loop */
    t5=clock();
    timespent=(float)(t5-t1)/CLOCKS_PER_SEC;
    sf_warning("CPU time for modeling and separation.: %f(second)",timespent);

    free(lmatxp);
    free(lmatzp);
    free(rmatxp);
    free(rmatzp);
    free(mmatxp);
    free(mmatzp);

    free(lmatxs);
    free(lmatzs);
    free(rmatxs);
    free(rmatzs);
    free(mmatxs);
    free(mmatzs);

    free(pp);
    free(qq);

    free(*p1);
    free(*p2);
    free(*p3);
    free(*q1);
    free(*q2);
    free(*q3);

    free(*vp);
    free(*vs);
    free(*ep);
    free(*de);
    free(*th);

    free(coeff_2dx);
    free(coeff_2dz);
    free(coeff_1dx);
    free(coeff_1dz);

    free(ijkx);
    free(ijkz);
    exit(0);
}
