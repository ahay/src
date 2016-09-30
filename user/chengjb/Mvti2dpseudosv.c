/* 2-D two-components wavefield modeling using pseudo-pure mode qSV-wave equation in VTI media.

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
#include "fdcoef.h"
#include "zero.h"
#include "kykxkztaper.h"
#include "kykxkz2yxz.h"
#include "clipsmthspec.h"

/* wave-mode separation operators */
#include "curlpoldevvtisv.h"
#include "filter2dsep.h"

/* wavefield propagators */
#include "fwpvtipseudosv.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
		int	ix, iz, ixf, izf, jx, jz, ixx, izz, i,j,im, jm,nx,nz,nxpad,nzpad,it,ii,jj;
		float   kxmax,kzmax;

        float   A, f0, t, t0, dx, dz, dxf, dzf, dt, dkx, dkz, dt2, div;
        int     mm, nvx, nvz, ns;
        int     hnkx, hnkz, nkx, nkz, nxz, nkxz;
        int     hnkx1, hnkz1, nkx1, nkz1;
        int     isx, isz, isxm, iszm; /*source location */

        int     nxf, nzf;
        int     nstep;            /* every nstep in spatial grids to calculate filters sparsely*/
        int     itaper;           /* tapering or not for spectrum of oprtator*/

        float   *coeff_1dx, *coeff_1dz, *coeff_2dx, *coeff_2dz; /* finite-difference coefficient */
        float   *kx, *kz, *kkx, *kkz, *kx2, *kz2, **taper;

        float **acx, **acz, **acxx, **aczz;        /* divergence operator of P-wave for a location */
        float **asx, **asz, **asxx, **aszz;        /* polarization operator of P-wave for a location */
        float **asvx, **asvz, **asvxx, **asvzz;    /* projection deviation operator of P-wave for a location */

        float ****ex, ****ez;                      /* operator for whole model for P-wave*/
        float **exx, **ezz;                        /* operator for constant model for P-wave*/

        float **vp0, **vs0, **epsi, **del;         /* velocity model */
        float **p1, **p2, **p3, **q1, **q2, **q3, **p3c, **q3c, **sum;  /* wavefield array */

        clock_t t1, t2, t3, t4, t5;
        float   timespent; 
        float   fx, fz;
        char    *tapertype;

        int     isep=1;
        int     ihomo=1;

		double  vp2, vs2, ep2, de2;

        sf_init(argc,argv);

        sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6, Fo7, Fo8;
        sf_file Fo9, Fo10, Fo11, Fo12, Fo13, Fo14, Fo15, Fo16;
       
        t1=clock();
 
        /*  wavelet parameter for source definition */
        f0=30.0;                  
        t0=0.04;                  
        A=1.0;                  

        /* time samping paramter */
        if (!sf_getint("ns",&ns)) ns=301;
        if (!sf_getfloat("dt",&dt)) dt=0.001;
        if (!sf_getint("isep",&isep)) isep=0;             /* if isep=1, separate wave-modes */ 
        if (!sf_getint("ihomo",&ihomo)) ihomo=0;          /* if ihomo=1, homogeneous medium */ 
        if (!sf_getint("itaper",&itaper)) itaper=0;          /* if itaper=1, taper the wavenumber domain p=operators*/ 
        if (NULL== (tapertype=sf_getstring("tapertype"))) tapertype="D"; /* taper type*/ 
        if (!sf_getint("nstep",&nstep)) nstep=2;

        sf_warning("isep=%d",isep);
        sf_warning("ihomo=%d",ihomo);
        sf_warning("itaper=%d",itaper);
        sf_warning("tapertype=%s",tapertype);
        sf_warning("nstep=%d",nstep);

        sf_warning("ns=%d dt=%f",ns,dt);
        sf_warning("read velocity model parameters");

        /* setup I/O files */
        sf_file Fvp0, Fvs0, Feps, Fdel;

        Fvp0 = sf_input ("in");  /* vp0 using standard input */
        Fvs0 = sf_input ("vs0");  /* vs0 */
        Feps = sf_input ("epsi");  /* epsi */
        Fdel = sf_input ("del");  /* delta */

        /* Read/Write axes */
        sf_axis az, ax;
        az = sf_iaxa(Fvp0,1); nvz = sf_n(az); dz = sf_d(az)*1000.0;
        ax = sf_iaxa(Fvp0,2); nvx = sf_n(ax); dx = sf_d(ax)*1000.0;
        fx=sf_o(ax)*1000.0;
        fz=sf_o(az)*1000.0;

        /* source definition */
        isx=nvx/2;
        isz=nvz/2;
        //isz=nvz*2/5;

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

        t2=clock();

        /* setup I/O files */
        Fo1 = sf_output("out"); /* pseudo-pure qSV-wave x-component */
        Fo2 = sf_output("PseudoPureSVz"); /* pseudo-pure qSV-wave z-component */
        Fo3 = sf_output("PseudoPureSV"); /* scalar qSV-wave field using divergence operator */

        puthead3(Fo1, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, dt*(ns-1));
        puthead3(Fo2, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, dt*(ns-1));
        puthead3(Fo3, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, dt*(ns-1));

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
           hnkx=800.0/dxf;
           hnkz=800.0/dzf;
           nkx=2*hnkx+1;   /* operator length in kx-direction */
           nkz=2*hnkz+1;   /* operator length in kz-direction */

           /* truncated spatial operators length for filtering*/
           hnkx1=400.0/dxf;
           hnkz1=400.0/dzf;
           nkx1=2*hnkx1+1;
           nkz1=2*hnkz1+1;

           sf_warning("nx=%d nz=%d nxf=%d nzf=%d", nx,nz,nxf,nzf);
           sf_warning("dx=%f dz=%f dxf=%f dzf=%f", dx,dz,dxf,dzf);
           sf_warning("hnkx=%d hnkz=%d nkx=%d nkz=%d", hnkx, hnkz, nkx, nkz);
           sf_warning("hnkx1=%d hnkz1=%d nkx1=%d nkz1=%d", hnkx1, hnkz1, nkx1, nkz1);

           dkx=2*SF_PI/dxf/nkx;
           dkz=2*SF_PI/dzf/nkz;
		   kxmax=SF_PI/dxf;
		   kzmax=SF_PI/dzf;

           kx=sf_floatalloc(nkx);
           kz=sf_floatalloc(nkx);
           kkx=sf_floatalloc(nkx);
           kkz=sf_floatalloc(nkx);
           kx2=sf_floatalloc(nkx);
           kz2=sf_floatalloc(nkx);

           taper=sf_floatalloc2(nkz, nkx);

           // define axis samples and taper in wavenumber domain 
           kxkztaper(kx, kz, kkx, kkz, kx2, kz2, taper, nkx, nkz, hnkx, hnkz, dkx, dkz, kxmax, kzmax, tapertype);

		   p3c=sf_floatalloc2(nz,nx);
		   q3c=sf_floatalloc2(nz,nx);
           sum=sf_floatalloc2(nz,nx);

           if(ihomo==1)
           {
			  exx=sf_floatalloc2(nkz1, nkx1);
			  ezz=sf_floatalloc2(nkz1, nkx1);
           }
           else{ // to store spatail varaied operators
			  ex=sf_floatalloc4(nkz1, nkx1, nzf, nxf);
			  ez=sf_floatalloc4(nkz1, nkx1, nzf, nxf);
           }
           nkxz=nkx*nkz;

           acx=sf_floatalloc2(nkz, nkx);
           acz=sf_floatalloc2(nkz, nkx);
           acxx=sf_floatalloc2(nkz, nkx);
           aczz=sf_floatalloc2(nkz, nkx);

           asx=sf_floatalloc2(nkz, nkx);
           asz=sf_floatalloc2(nkz, nkx);
           asxx=sf_floatalloc2(nkz, nkx);
           aszz=sf_floatalloc2(nkz, nkx);

		   asvx=sf_floatalloc2(nkz, nkx);
		   asvz=sf_floatalloc2(nkz, nkx);
		   asvxx=sf_floatalloc2(nkz, nkx);
		   asvzz=sf_floatalloc2(nkz, nkx);

           /* setup I/O files */
           Fo4 = sf_output("acx"); /* divergence x-comp */
           Fo5 = sf_output("acz"); /* divergence z-comp */
           Fo6 = sf_output("asx"); /*  qSV-wave's polarization x-comp */
           Fo7 = sf_output("asz"); /*  qSV-wave's polarization z-comp */
           Fo8 = sf_output("asvx"); /* qSV-wave projection deviation x-comp */
           Fo9 = sf_output("asvz"); /* qSV-wave projection deviation z-comp */

           Fo10 = sf_output("acxx"); /* curl operator x-comp in (x,z) domain */
           Fo11 = sf_output("aczz"); /* curl operator z-comp in (x,z) domain */
           Fo12 = sf_output("asxx");  /* qSV-wave's polarization x-comp in (x,z) domain */
           Fo13 = sf_output("aszz");  /* qSV-wave's polarization z-comp in (x,z) domain */
           Fo14 = sf_output("asvxx"); /* qSV-wave projection deviation x-comp in (x,z) domain */
           Fo15 = sf_output("asvzz"); /* qSV-wave projection deviation z-comp in (x,z) domain */

           puthead1(Fo4, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
           puthead1(Fo5, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
           puthead1(Fo6, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
           puthead1(Fo7, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
           puthead1(Fo8, nkz, nkx, dkz, -kzmax, dkx, -kxmax);
           puthead1(Fo9, nkz, nkx, dkz, -kzmax, dkx, -kxmax);

           puthead2(Fo10, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);
           puthead2(Fo11, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);
           puthead2(Fo12, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);
           puthead2(Fo13, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);
           puthead2(Fo14, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);
           puthead2(Fo15, nkz, nkx, dzf/1000.0, 0.0, dxf/1000.0, 0.0);

	   /*************calculate projection deviation grid-point by grid-point **********/
           for(ix=0,ixf=0;ix<nx;ix+=nstep,ixf++) 
           {
              for(iz=0,izf=0;iz<nz;iz+=nstep,izf++) 
              {
				vp2=vp0[ix][iz]*vp0[ix][iz];
				vs2=vs0[ix][iz]*vs0[ix][iz];
				ep2=1.0+2*epsi[ix][iz];
				de2=1.0+2*del[ix][iz];

                if(ixf%50==0&&izf%100==0) sf_warning("Operator: nxf=%d ixf=%d izf=%d vp2=%f vs2=%f",nxf, ixf,izf,vp2,vs2);

   	        /*************calculate projection deviation without tapering **********/
                //sf_warning("calculate projection deviation operators");
                /* devvtip: projection deviation operators for P-wave in VTI media */
                //divpoldevvtip(acx,acz,asx,asz,asvx,asvz,kx,kz,kkx,kkz,kx2,kz2,taper,hnkx,hnkz,
                              //dkx,dkz,vp2,vs2,ep2,de2,itaper);
							  
                curlpoldevvtisv(acx,acz,asx,asz,asvx,asvz,kx,kz,kkx,kkz,kx2,kz2,taper,hnkx,hnkz,
                              dkx,dkz,vp2,vs2,ep2,de2,itaper);

                /* inverse Fourier transform */
                kxkz2xz(asvx, asvxx, hnkx, hnkz, nkx, nkz);
                kxkz2xz(asvz, asvzz, hnkx, hnkz, nkx, nkz);

                // truncation and saving of operator in space-domain
                if(ihomo==1)
                {
                   for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++) 
                   for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++) 
                   {
                     exx[jx+hnkx1][jz+hnkz1]=asvxx[ixx][izz]; 
                     ezz[jx+hnkx1][jz+hnkz1]=asvzz[ixx][izz]; 
                   }
                }else{
                   for(jx=-hnkx1,ixx=hnkx-hnkx1;jx<=hnkx1;jx++,ixx++) 
                   for(jz=-hnkz1,izz=hnkz-hnkz1;jz<=hnkz1;jz++,izz++) 
                   {
                     ex[ixf][izf][jx+hnkx1][jz+hnkz1]=asvxx[ixx][izz]; 
                     ez[ixf][izf][jx+hnkx1][jz+hnkz1]=asvzz[ixx][izz]; 
                   }
                }
                if((ixf==nxf/2&&izf==nzf/2&&ihomo==0)||ihomo==1)
                {
                   //sf_warning("write-disk projection-deviation operators in kx-kz domain");
                   ikxkz2xz(acx, acxx, hnkx, hnkz, nkx, nkz);
                   ikxkz2xz(acz, aczz, hnkx, hnkz, nkx, nkz);
                   ikxkz2xz(asx, asxx, hnkx, hnkz, nkx, nkz);
                   ikxkz2xz(asz, aszz, hnkx, hnkz, nkx, nkz);

				   sf_floatwrite(acx[0], nkxz, Fo4);
				   sf_floatwrite(acz[0], nkxz, Fo5);
				   sf_floatwrite(asx[0], nkxz, Fo6);
				   sf_floatwrite(asz[0], nkxz, Fo7);
				   sf_floatwrite(asvx[0], nkxz, Fo8);
				   sf_floatwrite(asvz[0], nkxz, Fo9);

				   sf_floatwrite(acxx[0], nkxz, Fo10);
				   sf_floatwrite(aczz[0], nkxz, Fo11);
				   sf_floatwrite(asxx[0], nkxz, Fo12);
				   sf_floatwrite(aszz[0], nkxz, Fo13);
				   sf_floatwrite(asvxx[0], nkxz, Fo14);
				   sf_floatwrite(asvzz[0], nkxz, Fo15);
                }
                if(ihomo==1) goto loop;

              }// iz loop
            }//ix loop
            loop:;
            free(*acx);
            free(*acz);
            free(*acxx);
            free(*aczz);
            free(*asx);
            free(*asz);
            free(*asxx);
            free(*aszz);
            free(*asvx);
            free(*asvz);
            free(*asvxx);
            free(*asvzz);

            free(*taper);

            free(kx);
            free(kz);
            free(kx2);
            free(kz2);
            free(kkx);
            free(kkz);
       }//isep==1
       /****************End of Calculating Projection Deviation Operator****************/
       t3=clock();
       timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
       sf_warning("Computation time (operators): %f (second)",timespent);

       /****************begin to calculate wavefield****************/
       /****************begin to calculate wavefield****************/
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
        
        coeff1d(coeff_1dx,dx);
        coeff1d(coeff_1dz,dz);

        if(isep==1)
        {
             Fo16 = sf_output("PseudoPureSepSV"); /* scalar P-wave field using polarization projection oprtator*/

             puthead3(Fo16, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz/1000.0, fx/1000.0, dt*(ns-1));
        }

        sf_warning("==================================================");
        sf_warning("==  Porpagation Using Pseudo-Pure P-Wave Eq.    ==");
        sf_warning("==================================================");
        t4=clock();

	for(it=0;it<ns;it++)
	{
		t=it*dt;

		p2[isxm][iszm]+=Ricker(t, f0, t0, A);
		q2[isxm][iszm]+=Ricker(t, f0, t0, A);

                /* fwpvtipseudop: forward-propagating in VTI media with pseudo-pure P-wave equation */
		fwpvtipseudosv(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz,
                              nx, nz, vp0, vs0, epsi, del);

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

        t4=clock();
        timespent=(float)(t4-t3)/CLOCKS_PER_SEC;
        sf_warning("Computation time (propagation): %f (second)",timespent);
            if(isep==1)
            {
            //////////////////////////////////////////////////////////////////////////////////////////
            /* correct projection error for accurate separate qP wave in spatial domain */
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

			sf_floatwrite(sum[0],nx*nz, Fo16);
          
            t5=clock();
            timespent=(float)(t5-t4)/CLOCKS_PER_SEC;
            sf_warning("Computation time (one-step separation): %f (second)",timespent);
            } // isep==1
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

	exit(0);
}
