/* 2D finite difference modeling	*/

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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
	int it, ix, iz;			// index
	int sx, sz;			// shot positions
	int rx0, rz, nr, dr;			// recevier positions
	int nt, nx, nz, nwv;	 	// dimensions
	float dt, dx, dz, dt2;		// increments
	int wvlt, owv;				// arguments
	float wvp[4], *pwv=NULL;
	int jt;

	float **vv, **pout;
	float **u0,**u1, **ud, u2; 	// tmp arrays 
	void * h;
	sf_file data, modl, wave; /* I/O files */

	sf_init(argc,argv);

	modl = sf_input ("in");   //- velocity model 
	data = sf_output("out");  //- seismic data
 
	if(sf_getstring("wave")!=NULL) 
		wave = sf_output("wave"); /* wavefield movie file */
	else wave=NULL;

//	div  = sf_output("div"); /* wavefield movie file */
	
	if (!sf_getint("nt",&nt)) nt=1000; 
	/* time samples */
	if (!sf_getfloat("dt",&dt)) dt=0.004; 
	/* time interval */
	if (!sf_getint("jt",&jt)) jt=40; 
	/* wave movie time interval */


	nwv=0; owv=0;
	if (!sf_getint("wvlt", &wvlt)) wvlt=0; 
	/* wavelet type "ricker/harmonic/other" */
	if (!sf_getfloat("w0",&u2)) u2=35.0; 
	/* central frequency for ricker/harmonic wavelet */
	wvp[0] = u2;
	wvp[1] = 0.0;	// zero phase
	switch(wvlt)
	{
	case 0:		// ricker wavelet
//		owv = 2.0 /(u2*dt);
//		owv = (owv<nt)? -owv:-nt;
		nwv = nt- owv;
		pwv = sf_floatalloc(nwv);
		for(it=owv; it<nt; it++) pwv[it-owv] = dt*it;
		sf_wvlt_rck(nwv, pwv, wvp);
		break;
	case 1:		// sin signal
		nwv = nt;
		owv = 0;
		pwv = sf_floatalloc(nwv);
		for(it=owv; it<nt; it++) pwv[it-owv] = dt*it;
		sf_wvlt_harmonic(nwv, pwv, wvp);
		break;
	default:
		sf_error("wvlt = %d not support", wvlt);
	}

	if (!sf_histint(modl,"n1", &nz)) sf_error("n1");
	if (!sf_histint(modl,"n2", &nx)) sf_error("n2");
	if (!sf_histfloat(modl,"d1", &dz)) sf_error("d1");
	if (!sf_histfloat(modl,"d2", &dx)) sf_error("d2");

	if (!sf_getint("sx",&sx)) sx=0; 
	/* x position index of the source */
	if (!sf_getint("sz",&sz)) sz=0; 
	/* z position index of the source */
	if(sx < 0 || sx >= nx) 
		sf_error("sx (= %d) should between (0,%d)",sx,nx-1);
	if(sz < 0 || sz >= nz) 
		sf_error("sz (= %d) should between (0,%d)",sz,nz-1);


	if (!sf_getint("rx0",&rx0)) rx0=0; 
	/* x position index of first receiver */
	if (!sf_getint("nr",&nr)) nr=1; 
	/* receiver numbers */
	if (!sf_getint("dr",&dr)) dr=1; 
	/* receiver interval of unit "dx" */
	if (!sf_getint("rz",&rz)) rz=0; 
	/* z position index of receivers */
	if(rx0 < 0 || rx0+dr*(nr-1) >= nx) 
		sf_error("rx0+dr*(nr-1) should between (0,%d)",nx-1);
	if(rz < 0 || rz >= nz) 
		sf_error("rz (= %d) should between (0,%d)",rz,nz-1);
	

	sf_putint(data, "n1", nt);
	sf_putint(data, "n2", nr);
	sf_putfloat(data, "d1", dt);
	sf_putfloat(data, "d2", dr*dx);
	sf_putfloat(data, "o2", rx0*dx);

	if(wave!=NULL)
	{
		sf_putint(wave, "n1", nz);
		sf_putint(wave, "n2", nx);
		sf_putint(wave, "n3", nt/jt);
		sf_putfloat(wave, "d1", dz);
		sf_putfloat(wave, "d2", dx);
		sf_putfloat(wave, "d3", dt*jt);
		sf_putfloat(wave, "o3", 0.0);
	}

//	sf_putint(div, "n1", nz);
//	sf_putint(div, "n2", nx);
//	sf_putint(div, "n3", nwv);
//	sf_putfloat(div, "d1", dz);
//	sf_putfloat(div, "d2", dx);
//	sf_putfloat(div, "d3", dt);


	dt2 = dt*dt;
	
	h = sf_fd4_init(nz, nx, dz, dx);

	/* allocate temporary arrays */
	vv = sf_floatalloc2(nz, nx);
	pout = sf_floatalloc2(nt, nr);
	u0=sf_floatalloc2(nz,nx);		// previous data
	u1 = sf_floatalloc2(nz, nx);
	ud=sf_floatalloc2(nz,nx);	// updata

	sf_floatread(vv[0], nx*nz, modl);

	for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
			vv[ix][iz] = vv[ix][iz] * vv[ix][iz]*dt2;
			ud[ix][iz] = 0.0;
			u0[ix][iz] = 0.0;
			u1[ix][iz] = 0.0;
		}
	}

	/* Time loop */
	for (it=owv; it < nt; it++) 
	{
		sf_fd4_laplacian(h, u1, ud);
		ud[sx][sz] += pwv[it-owv];
		//		if(div!=NULL) sf_floatwrite(ud[0], nz*nx, div);
	
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
		}
		if(it>=0)
		{
			for(ix=0; ix<nr; ix++)
				pout[ix][it] = u0[ix*dr+rx0][rz];
		}
		if(wave!=NULL && it>=0 && it%jt == 0)// wave
		{
			sf_floatwrite(u1[0], nz*nx, wave);
		}
	}
	/* output seismic data */
	sf_floatwrite(pout[0], nt*nr, data);

	sf_fd4_close(h);
	
	free(pwv);
	free(vv[0]);
	free(u0[0]);
	free(u1[0]);
	free(pout[0]);
	free(ud[0]);
	free(vv);
	free(u0);
	free(u1);
	free(pout);
	free(ud);
	return (0);
}
