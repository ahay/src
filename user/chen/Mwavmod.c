/* 1-2-3D finite difference modeling	*/

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
#include "wavmod.h"
#include "fd3.h"

int main(int argc, char* argv[])
{
	int ix, iy, it;			// index
	int nz, nx, ny, nt, n1;		 	// dimensions
	float dz, dx, dy, dt;		// increments
	float oz, ox, oy, ot, o1;		// original

	sf_file hwvlt, hwave, hdata, hvel; /* I/O files */
	float *wvlt, ***wave, ***data;

	int xs, ys, zs;			// shot positions
	int zr;		// recevier positions
	int jt;

	sf_init(argc,argv);

	hwvlt = sf_input ("in");   //- source wavelet
	hvel  = sf_input ("vel");  //- velocity field
	hdata = sf_output("out");  //- seismic data
 
	if(sf_getstring("wave")!=NULL) 
		hwave = sf_output("wave"); /* wavefield movie file */
	else hwave=NULL;

	if (!sf_histint(hwvlt,"n1", &n1)) sf_error("n1 needed in wvlt");
	if (!sf_histfloat(hwvlt,"d1", &dt)) sf_error("d1 needed in wvlt");
	if (!sf_histfloat(hwvlt,"o1", &o1)) o1 = 0;

	if (!sf_histint(hvel,"n1", &nz)) sf_error("n1 needed in vel");
	if (!sf_histint(hvel,"n2", &nx)) nx = 1;
	if (!sf_histint(hvel,"n3", &ny)) ny = 1;
	if (!sf_histfloat(hvel,"d1", &dz)) sf_error("d1 needed in vel");
	if (!sf_histfloat(hvel,"d2", &dx)) dx = 0.0;
	if (!sf_histfloat(hvel,"d3", &dy)) dy = 0.0;
	if (!sf_histfloat(hvel,"o1", &oz)) oz = 0.0;
	if (!sf_histfloat(hvel,"o2", &ox)) ox = 0.0;
	if (!sf_histfloat(hvel,"o3", &oy)) oy = 0.0;

	if (!sf_getint("jt",&jt)) jt=100; 
	/* wave movie time interval */
	if (!sf_getfloat("ot", &ot)) ot = o1; 
	/* time delay */

	if (!sf_getint("xs",&xs)) xs=0; 
	if (!sf_getint("ys",&ys)) ys=0; 
	if (!sf_getint("zs",&zs)) zs=0; 
	/* x-y-z position index of the source */
	if(xs < 0 || xs >= nx) 
		sf_error("xs (= %d) should between (0,%d)",xs,nx-1);
	if(ys < 0 || ys >= ny) 
		sf_error("ys (= %d) should between (0,%d)",ys,ny-1);
	if(zs < 0 || zs >= nz) 
		sf_error("zs (= %d) should between (0,%d)",zs,nz-1);

	if (!sf_getint("zr",&zr)) zr=1; 
	/* z dimension of receiver */
	if(zr < 0 || zr >= nz) 
		sf_error("z position should between (0,%d)",nz-1);
	
	nt = (dt*n1+o1-ot)/dt +1;

	sf_putint(hdata, "n1", nt);
	sf_putfloat(hdata, "o1", ot);
	sf_putfloat(hdata, "d1", dt);
	sf_putint(hdata, "n2", nx);
	sf_putfloat(hdata, "o2", ox);
	sf_putfloat(hdata, "d2", dx);
	sf_putint(hdata, "n3", ny);
	sf_putfloat(hdata, "o3", oy);
	sf_putfloat(hdata, "d3", dy);

	if(wave!=NULL)
	{
		sf_putint(hwave, "n1", nz);
		sf_putint(hwave, "n2", nx);
		sf_putint(hwave, "n3", ny);
		sf_putint(hwave, "n4", nt/jt);
		sf_putfloat(hwave, "d1", dz);
		sf_putfloat(hwave, "d2", dx);
		sf_putfloat(hwave, "d3", dy);
		sf_putfloat(hwave, "d4", dt*jt);
		sf_putfloat(hwave, "o1", oz);
		sf_putfloat(hwave, "o2", ox);
		sf_putfloat(hwave, "o3", oy);
		sf_putfloat(hwave, "o4", o1);
	}

	fd3_init(dz, dx, dy);
	wavmod_init(nz, nx, ny, dt, hvel, zs, xs, ys);

	/* allocate temporary arrays */
	wvlt = sf_floatalloc(n1);
	sf_floatread(wvlt, n1, hwvlt);

	data = sf_floatalloc3(nt, nx, ny);
	wave  = sf_floatalloc3(nz, nx, ny);

	memset(wave[0][0], 0, nx*ny*nz);

	for (it=0; it < nt; it++) 
	{
		wavmod(wave[0][0], wvlt[it]);

		if(dt*it+o1 >= ot)
		{
			for(iy=0; iy<ny; iy++)
			for(ix=0; ix<nx; ix++)
				data[iy][ix][it] = wave[iy][ix][zr];
		}
		if(wave!=NULL && it>=0 && it%jt == 0)// wave
		{
			sf_floatwrite(wave[0][0], nz*nx*ny, hwave);
		}
	}
	/* output seismic data */
	sf_floatwrite(data[0][0], nt*ny*nx, hdata);

	wavmod_close();

	free(wave[0][0]);
	free(wave[0]);
	free(wave);
	free(data[0][0]);
	free(data[0]);
	free(data);
	free(wvlt);
	return (0);
}
