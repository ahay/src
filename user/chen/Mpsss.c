/* phase shift wave extrapolation. */

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

#include "psss.h"
#include "rfft2.h"


int main(int argc, char* argv[])
{
	int ix, iz, jz, njz;
	float dt,dx,dz;
	int nx,nz,nt,nz0,nw,nt1,nx1,rz;
	float ox, **ptx, **pim, *vel;  
	sf_complex **pfk;
	void *h;

	sf_file modl,data,wave,imag;

	sf_init(argc,argv);

	data = sf_input("in");  	/* seismic data, zero offset u(t,0,x) */
	imag = sf_output ("out"); 	/* image data: u(0,z,x) */
	modl = sf_input("vel"); 	/* velocity model m(z,x) */
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
	sf_putfloat(imag,"d1",dz*jz);
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

	nt1 = nt;
	nx1 = nx;
	h = sf_rfft2_init(&nt1, &nx1, &nw, 1);

	/* read data and velocity */
	vel = sf_floatalloc(nz0);
	sf_floatread(vel,nz0,modl);

	ptx  = sf_floatalloc2(nt,nx);	// U_z(t,x)
	pim  = sf_floatalloc2(nz,nx);		// u(x,z,0)
	pfk  = sf_complexalloc2(nw,nx1);

	sf_floatread(ptx[0],nt*nx,data);

	for(ix=0;ix<nx;ix++)	
		pim[ix][0]=ptx[ix][0];  // imag iz =0

	sf_rfft2(h, ptx, pfk);

	sf_floatwrite(ptx[0],nt*nx,wave); // wave slice iz=0

	sf_psss_init(nw,nx1,nz0,
		1.0/(dt*nt1),dx,dz,vel);


	for(iz=1;iz<nz;iz++)
	{
		sf_psss_step(iz-1+rz, pfk);

		sf_rifft2(h, pfk, ptx);
		for(ix=0;ix<nx;ix++)	
			pim[ix][iz]=ptx[ix][0];  // imag iz =0

		if((iz) % jz == 0)
			sf_floatwrite(ptx[0],nt*nx,wave); 
		sf_warning("%d;",iz);
	}

	sf_floatwrite(pim[0],nz*nx,imag);
 
	sf_psss_close();
	sf_rfft2_close(h);

	free(vel);
	free(ptx[0]);
	free(pim[0]);
	free(pfk[0]);
	free(ptx);
	free(pim);
	free(pfk);
	return 0;
}

