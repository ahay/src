/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in VTI media.
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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
#include "_cjb.h"

/* calculate phase, group velocity and angle  */
#include "zero.h"
#include "ricker.h"
#include "puthead.h"
#include "phasegroup2dvti.h"

int main(int argc, char* argv[])
{

        int    nx, nz;
        float  dx, dz, time, da; 
        float  vp0, vs0, eps, del, the;
        float  f0, t0;

        sf_init(argc,argv);

        if (!sf_getint("nx",&nx)) nx=201;
        if (!sf_getint("nz",&nz)) nz=201;
        if (!sf_getfloat("dx",&dx)) dx=0.008;
        if (!sf_getfloat("dz",&dz)) dz=0.008;
        if (!sf_getfloat("time",&time)) time=0.2;  //unit: SECOND
        if (!sf_getfloat("da",&da)) da=0.05;
        if (!sf_getfloat("vp0",&vp0)) vp0=3000.0;
        if (!sf_getfloat("vs0",&vs0)) vs0=1200.0;
        if (!sf_getfloat("eps",&eps)) eps=0.2;
        if (!sf_getfloat("del",&del)) del=0.1;
        if (!sf_getfloat("the",&the)) the=0.0;
        if (!sf_getfloat("t0",&t0)) t0=0.04;
        if (!sf_getfloat("f0",&f0)) f0=20.0;

        sf_warning("nx= %d nz= %d",nx,nz);
        sf_warning("dx= %f dz= %f da= %f",dx,dz,da);
        sf_warning("time= %f ",time);
        sf_warning("vp0= %f vs0= %f eps= %f del=%f the=%f",vp0,vs0,eps,del,the);
        sf_warning("t0= %f f0= %f",t0,f0);

        int   it, ns=100;
        float x, a, dt=time/ns;
        
        float amax=0.0;
        int   itmax=0;
        for(it=0;it<ns;it++)
        {
          x=pow(SF_PI*f0*(it*dt-t0),2);
          a=-exp(-x)*(1-2*x);
          if(fabs(a)>amax){
             amax=fabs(a);
             itmax=it;
          } 
        }
        time -= itmax*dt;
        sf_warning("time= %f ",time);

	float **wfp, **wfs, **wf;
  
        wfp = sf_floatalloc2(nz,nx);
        wfs = sf_floatalloc2(nz,nx);
        wf = sf_floatalloc2(nz,nx);

	zero2float(wfp,nz,nx);
	zero2float(wfs,nz,nx);
	zero2float(wf,nz,nx);

	int   i,j,d_x,d_z;
	float vpp, vpg, apg, vsp, vsg, asg, ap, d;

        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++){
          wf[i][j] = 0.0;
          wfp[i][j] = 1;
          wfs[i][j] = -1;
        }

        /* setup I/O files */
        sf_file Fo1, Fo2, Fo3;
        Fo1 = sf_output("out"); // P & SV waves' wavefront
        Fo2 = sf_output("WFp"); // P-wave's wavefront
        Fo3 = sf_output("WFs"); // SV-wave's wavefront

        puthead2(Fo1, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo2, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo3, nz, nx, dz, 0.0, dx, 0.0);

        // source_x source_z
        int sx=nx/2;
        int sz=nz/2; 

        dx *=1000.0;
        dz *=1000.0;

	for(i=0;i<360/da;i++)
	{
		ap = i*da*SF_PI/180.0;
		vpp = vpphase2dvti(vp0, vs0, eps, del, ap);
		vsp = vsphase2dvti(vp0, vs0, eps, del, ap);

		vapgroup2dvti(vp0, vs0, eps, del, ap, &vpg, &apg);
		vasgroup2dvti(vp0, vs0, eps, del, ap, &vsg, &asg);

		if(i>=180/da){
		   apg+=SF_PI;
		   asg+=SF_PI;
                }
                //sf_warning("vpp=%f vsp=%f vpg=%f vsg=%f apg=%f asg=%f",vpp,vsp,vpg,vsg,apg,asg);

		d=vpg*time;
		d_x=(int)(d*sin(the+apg)/dx)+sx;
		d_z=-(int)(d*cos(the+apg)/dz)+sz;
	
                //sf_warning("P:  d_x=%d d_z=%d",d_x,d_z);

		if(d_x>=0 && d_x<nx && d_z>=0 && d_z<nz)
			wfp[d_x][d_z]=-1;

		d=vsg*time;
		d_x=(int)(d*sin(the+asg)/dx)+sx;
		d_z=-(int)(d*cos(the+asg)/dz)+sz;
	
                //sf_warning("S:  d_x=%d d_z=%d",d_x,d_z);
		if(d_x>=0 && d_x<nx && d_z>=0 && d_z<nz)
			wfs[d_x][d_z]=1;
        }

        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++){
           wf[i][j] = wfp[i][j]+wfs[i][j];
           sf_floatwrite(&wf[i][j], 1, Fo1);
           sf_floatwrite(&wfp[i][j], 1, Fo2);
           sf_floatwrite(&wfs[i][j], 1, Fo3);
        }
/*
        FILE *fp;
        fp=fopen("WF","wb");
        fwrite(&wf[0][0], sizeof(float), nx*nz, fp);
        fclose(fp);
*/
        free(*wfp);
        free(*wfs);
        free(*wf);
}
