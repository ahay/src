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

#include "psefd.h"
#include "rfft1.h"


int main(int argc, char* argv[])
{
	int it,iw,ix,iz,jz,njz;
	float dt,dx,dz;
	int nw,nx,nz,nt,nz0,nt2,rz;

	float *v1; 
	sf_complex **u1; 

	float ox;
	float **u2,**u3,**vel;  
	void * h;

	sf_file modl,data,wave,imag;

	sf_init(argc,argv);

    data = sf_input("in");  	/* seismic data, zero offset u(t,0,x) */
    imag = sf_output ("out"); 	/* image data: u(0,z,x) */
	modl = sf_input("vel"); 		/* velocity model m(z,x) */
    wave = sf_output("wave"); 	/* wavefield: u(t,z,x) */

    if (!sf_histint(modl,"n1",&nz0)) sf_error("n1");
    if (!sf_histint(modl,"n2",&nx)) sf_error("n2");
    if (!sf_histfloat(modl,"o2",&ox)) sf_error("o2");
    if (!sf_histint(data,"n1",&nt)) sf_error("n1");
    if (!sf_histfloat(data,"d1",&dt)) sf_error("d2");

    if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
    if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");

    if (!sf_getint("jz",&jz)) jz=40; 
    /* depth step for wave data */
    if (!sf_getint("rz",&rz)) rz=0; 
    /* receiver depth */
    if (!sf_getint("nz",&nz)) nz=nz0-rz; 
    /* depth number */

    sf_putint(imag,"n1",nz);
    sf_putint(imag,"n2",nx);
    sf_putfloat(imag,"d1",dz);
    sf_putfloat(imag,"d2",dx);
    sf_putfloat(imag,"o1",dz*rz);
    sf_putfloat(imag,"o2",ox);

	njz=1+(nz-1)/jz;
    sf_putint(wave,"n1",nt);
    sf_putint(wave,"n2",nx);
    sf_putint(wave,"n3",njz);
    sf_putfloat(wave,"d1",dt);
    sf_putfloat(wave,"d2",dx);
    sf_putfloat(wave,"d3",dz*jz);
    sf_putfloat(wave,"o1",0);
    sf_putfloat(wave,"o2",ox);
    sf_putfloat(wave,"o3",dz*rz);

	nt2 = nt;
	h = sf_rfft1_init(&nt2, &nw, 1);
    v1 = sf_floatalloc(nt);	

    /* read data and velocity */
    vel = sf_floatalloc2(nz0,nx);
    sf_floatread(vel[0],nz0*nx,modl);

    u1 = sf_complexalloc2(nw,nx);	// U_z(w,x)
    u2 = sf_floatalloc2(nt,nx);	// u_z(x,t)
    u3 = sf_floatalloc2(nz,nx);		// u(x,z,0)


	for(ix=0;ix<nx;ix++)	
	{
    	sf_floatread(v1,nt,data);
		for(it=0;it<nt;it++)   u2[ix][it]=v1[it];
        u3[ix][0]=v1[0];  // imag rz

		sf_rfft1(h, v1, u1[ix]);
	}
	sf_floatwrite(u2[0],nt*nx,wave); // wave slice iz=0

	sf_psefd_init(nx,nw,dx,dz,
		1.0/(nt2*dt),vel);


	for(iz=1;iz<nz;iz++)
	{
		sf_psefd_step3(iz+rz,u1);

		if(iz % jz == 0 )
		{
			for(ix=0;ix<nx;ix++)
			{
				sf_rifft1(h, u1[ix], v1);
				for(it=0;it<nt;it++)	u2[ix][it]=v1[it];
				u3[ix][iz]=v1[0];
			}
			sf_floatwrite(u2[0],nt*nx,wave); 
		} else {
			for(ix=0;ix<nx;ix++)
			{
				u3[ix][iz]=0.0;
				for(iw=0;iw<nw-1;iw++)
					u3[ix][iz] += creal(u1[ix][iw]);
				u3[ix][iz]*=2.0/nt2;
			}
		}

		sf_warning("%d;",iz);
	}

	sf_floatwrite(u3[0],nz*nx,imag);
 
	sf_psefd_close();

	return 0;
}

