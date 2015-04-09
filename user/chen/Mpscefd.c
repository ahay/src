/* EFD phase shift wave extrapolation. */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

#include "pscefd.h"


int main(int argc, char* argv[])
{
	int iw,it,ix,iz,jt,njt;
	float dt,dx,dz;
	int nw,nx,nz,nt,nz0,nt2;

	kiss_fftr_cfg cfg,icfg;
	float *v1; 		/* fft buffer */
	kiss_fft_cpx *v2,**u1; /* fft buffer */

	float ox;
	float **u2,**u3,**vel;  

	sf_file modl,data,wave,imag;

	sf_init(argc,argv);

	modl = sf_input("in"); 		/* velocity model m(z,x) */
    imag = sf_output ("out"); 	/* image data: u(0,z,x) */
    data = sf_input("data");  	/* seismic data, zero offset u(t,0,x) */
    wave = sf_output("wave"); 	/* wavefield: u(t,z,x) */

    if (!sf_histint(modl,"n1",&nz0)) sf_error("n1");
    if (!sf_histint(modl,"n2",&nx)) sf_error("n2");
    if (!sf_histfloat(modl,"o2",&ox)) sf_error("o2");
    if (!sf_histint(data,"n1",&nt)) sf_error("n1");
    if (!sf_histfloat(data,"d1",&dt)) sf_error("d2");

    if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
    if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");

    if (!sf_getint("jt",&jt)) jt=40; 
    /* time step for wave data */
    if (!sf_getint("nz",&nz)) nz=nz0; 
    /* depth number */

    sf_putint(imag,"n1",nz);
    sf_putint(imag,"n2",nx);
    sf_putfloat(imag,"o2",ox);

	njt=1+(nt-1)/jt;
    sf_putint(wave,"n1",njt);
    sf_putfloat(wave,"d1",jt*dt);
    sf_putfloat(wave,"o1",0);
    sf_putint(wave,"n2",nx);
    sf_putfloat(wave,"d2",dx);
    sf_putfloat(wave,"o2",ox);
    sf_putint(wave,"n3",nz);
    sf_putfloat(wave,"d3",dz);
    sf_putfloat(wave,"o3",0);

	
	nw=kiss_fft_next_fast_size((nt+1)/2);
    nt2 = 2*(nw-1);

    cfg  = kiss_fftr_alloc(nt2,0,NULL,NULL);
    icfg = kiss_fftr_alloc(nt2,1,NULL,NULL);

    /* read data and velocity */
    vel = sf_floatalloc2(nz,nx);
    sf_floatread(vel[0],nz*nx,modl);

    v1 = sf_floatalloc(nt2);
    v2 = (kiss_fft_cpx *) sf_complexalloc(nw);

    u1 = (kiss_fft_cpx **) sf_complexalloc2(nx,nw);	/* U_z(x,w) */
    u2 = sf_floatalloc2(njt,nx);	/* u_z(x,t) */
    u3 = sf_floatalloc2(nz,nx);		/* u(x,z,0) */

	for(it=nt;it<nt2;it++) v1[it]=0.0;

	for(ix=0;ix<nx;ix++)	
	{
    	sf_floatread(v1,nt,data);
		for(it=0;it<njt;it++)   u2[ix][it]=v1[it*jt];
		u3[ix][0]=v1[0];  /* imag iz =0 */

		kiss_fftr(cfg,v1,v2);
		kiss_fftri(icfg,v2,v1);
		for(iw=0;iw<nw;iw++) 
		{
			u1[iw][ix].r=v2[iw].r; 
			u1[iw][ix].i=v2[iw].i; 
		}
	}
	sf_floatwrite(u2[0],njt*nx,wave); /* wave slice iz=0 */

/*	sf_pscefd_init(nx,nw,dx,dz,
//		1.0/(nt2*dt),vel); */


	for(iz=1;iz<nz;iz++)
	{
/*		sf_pscefd_step3(iz,u1); */

		for(ix=0;ix<nx;ix++)
		{
			for(iw=0;iw<nw;iw++)
			{
				v2[iw].r=u1[iw][ix].r;
				v2[iw].i=u1[iw][ix].i;
			}
			kiss_fftri(icfg,v2,v1);
			for(it=0;it<njt;it++)	u2[ix][it]=v1[it*jt]/nt2;
			u3[ix][iz]=v1[0]/nt2;
		}
		sf_floatwrite(u2[0],njt*nx,wave); 

		sf_warning("%d;",iz);
	}

	sf_floatwrite(u3[0],nz*nx,imag);
 
/*	sf_pscefd_exit(); */

	return 0;
}

