/* 2D finite difference modeling	*/
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
#include "fd4.h"
#include "wavelet.h"

int main(int argc, char* argv[])
{
	int it, ix, iz;		// index
	int sx, sy, sz;		// shot positions
	int nt, nx, nz, nwv;	 		// dimensions
	float dt, dx, dz, dt2;				// increments
	int wvlt, owv;						// arguments
	float wvp[4], *pwv;
	int jt, rz;

	float **vv, **pout;
	float **u0,**u1, **ud, u2; 	// tmp arrays 
	void * h;
	sf_file data, modl, wave, div; /* I/O files */

	sf_init(argc,argv);

	modl = sf_input ("in");   //- velocity model 
	data = sf_output("out");  //- seismic data
 
	wave = sf_output("wave"); /* wavefield movie file */
	div  = sf_output("div"); /* wavefield movie file */
	
	if (!sf_getint("nt",&nt)) nt=1000; 
	/* time samples */
	if (!sf_getfloat("dt",&dt)) dt=0.004; 
	/* time interval */
	if (!sf_getint("jt",&jt)) jt=20; 
	/* wave movie time interval */
	if (!sf_getint("rz",&rz)) rz=0; 
	/* reciever depth */

	if (!sf_getint("sx",&sx)) sx=0; 
	/* shot x position for 2/3D modeling  */
	if (!sf_getint("sy",&sy)) sy=0; 
	/* shot y position for 3D modeling */
	if (!sf_getint("sz",&sz)) sz=0; 
	/* shot z position for 2/3D modeling */

	nwv=0; owv=0;
	if (!sf_getint("wvlt", &wvlt)) wvlt=0; 
	/* wavelet type "ricker/other" */
	switch(wvlt)
	{
	case 0:
		if (!sf_getfloat("w0",&u2)) u2=35.0; 
		/* central frequency of ricker wavelet */
		wvp[0] = u2;
		wvp[1] = 0.0;
		owv = 5.0 /(u2*dt);
		owv = (owv<nt)? -owv:-nt;
		nwv = nt- owv;
		break;
	case 1:
		sf_error("wvlt = %d not support", wvlt);
	}

	if (!sf_histint(modl,"n1", &nz)) sf_error("n1");
	if (!sf_histint(modl,"n2", &nx)) sf_error("n2");
	if (!sf_histfloat(modl,"d1", &dz)) sf_error("d1");
	if (!sf_histfloat(modl,"d2", &dx)) sf_error("d2");


	pwv = sf_floatalloc(nwv);
	for(it=owv; it<nt; it++) pwv[it-owv] = dt*it;
	sf_wvlt_rck(nwv, pwv, wvp);

	sf_putint(data, "n1", nt);
	sf_putint(data, "n2", nx);
	sf_putfloat(data, "d1", dt);
	sf_putfloat(data, "d2", dx);

	sf_putint(wave, "n1", nz);
	sf_putint(wave, "n2", nx);
	sf_putint(wave, "n3", nt/jt);
	sf_putfloat(wave, "d1", dz);
	sf_putfloat(wave, "d2", dx);
	sf_putfloat(wave, "d3", dt*jt);


	dt2 = dt*dt;
	
	h = sf_fd4_init(nz, nx, dz, dx);

	/* allocate temporary arrays */
	vv = sf_floatalloc2(nz, nx);
	pout = sf_floatalloc2(nt, nx);
	u0=sf_floatalloc2(nz,nx);		// previous data
	u1 = sf_floatalloc2(nz, nx);
	ud=sf_floatalloc2(nz,nx);	// updata

	sf_floatread(vv[0], nx*nz, modl);

	for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
			ud[ix][iz] = 0.0;
			u0[ix][iz] = 0.0;
			u1[ix][iz] = 0.0;
			vv[ix][iz] = vv[ix][iz] * vv[ix][iz]*dt2;
		}
	}

	/* Time loop */
	for (it=owv; it < nt; it++) 
	{
		sf_fd4_laplacian(h, u1, ud);
		ud[sx][sz] += pwv[it-owv];
		if(div!=NULL) sf_floatwrite(ud[0], nz*nx, div);

		for (ix=0; ix<nx; ix++) 
		{
			for (iz=0; iz<nz; iz++) 
			{
				/* scale by velocity */
				ud[ix][iz] *= vv[ix][iz];
				/* time step */
				u2 = 2.0*u1[ix][iz] - u0[ix][iz] + ud[ix][iz]; 
		
				u0[ix][iz] = u1[ix][iz];
				u1[ix][iz] = u2;
			}
			if(it>=0) pout[ix][it] = u0[ix][rz];
		}
		if(it>=0 && it%jt == 0)// wave
		{
			sf_floatwrite(u1[0], nz*nx, wave);
		}
	}

	/* output seismic data */
	sf_floatwrite(pout[0], nt*nx, data);

	sf_fd4_release(h);
	
	free(pwv);
//	sf_free(vv[0]);
//	sf_free(u0[0]);
//	sf_free(u1[0]);
//	sf_free(pout[0]);
//	sf_free(ud[0]);
	return (0);
}
