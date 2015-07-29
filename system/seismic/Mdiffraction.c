/* Generate diffractions in zero-offset data. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "ricker.h"

int main(int argc, char *argv[])
{
    int nx, ny, nt, four, nsp, nd, ix, iy, isp, kz, kx, ky, k;
    float dx, dy, dt, t0, z, x, y, x1, y1, x0, y0, tx, ty;
    float freq, s1, s2, s12, t;
    float *v1, *v2, *v12, **sp, *trace, *time, *delt, *ampl;
    sf_file w1, w2, w12, spikes, data;

    sf_init(argc,argv);
    w1 = sf_input("in");
    w2 = sf_input("w2");
    w12 = sf_input("w12");
    spikes = sf_input("spikes");
    data = sf_output("out");
    
    if (!sf_histint(w1,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(w1,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(w1,"o1",&t0)) t0=0.0;

    if (!sf_histint(w1,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(w1,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(w1,"o2",&x0)) x0=0.0;

    if (!sf_histint(w1,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(w1,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(w1,"o3",&y0)) y0=0.0;

    if (!sf_histint(spikes,"n1",&four) || four != 4) sf_error("Need n1=4 in spikes");
    if (!sf_histint(spikes,"n2",&nsp)) nsp=1;

    nd = nt*nx*ny;

    v1 = sf_floatalloc(nd);
    v2 = sf_floatalloc(nd);
    v12 = sf_floatalloc(nd);
    
    sf_floatread(v1,nd,w1);
    sf_floatread(v2,nd,w2);
    sf_floatread(v12,nd,w12);

    sp = sf_floatalloc2(4,nsp);
    trace = sf_floatalloc(nt);

    sf_floatread(sp[0],4*nsp,spikes);

    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    time = sf_floatalloc(nsp);
    delt = sf_floatalloc(nsp);
    ampl = sf_floatalloc(nsp);

    sf_aastretch_init (false, nt, t0, dt, nsp);

    for (iy=0; iy < ny; iy++) {
	y1 = y0 + iy*dy;
	for (ix=0; ix < nx; ix++) {
	    x1 = x0 + ix*dx;
	    for (isp=0; isp < nsp; isp++) {
		z = sp[isp][0];
		x = sp[isp][1];
		y = sp[isp][2];

		kz = 0.5+(z-t0)/dt;
		kx = 0.5+(x-x0)/dx;
		ky = 0.5+(y-y0)/dy;

		x -= x1;
		y -= y1;

		if (kz < 0 || kz >= nt ||
		    kx < 0 || kx >= nx ||
		    ky < 0 || ky >= ny) {
		    time[isp] = t0-10*dt;
		    delt[isp] = 0.;
		    ampl[isp] = 0.0;
		} else {
		    k = kz+nt*(kx+nx*ky);

		    s1 = v1[k];
		    s2 = v2[k];
		    s12 = v12[k];
		    
		    t = sqrtf(z*z+4.*(s1*x*x+s2*y*y+2.*s12*x*y));
		    tx = 4.*(s1*x+s12*y)/(t+dt);
		    ty = 4.*(s2*y+s12*x)/(t+dt);

		    time[isp] = t;
		    delt[isp] = SF_MAX(fabsf(tx*dx),fabsf(ty*dy));
		    ampl[isp] = sp[isp][3];
		}
	    }

	    sf_aastretch_define (time,delt,NULL);
	    sf_aastretch_lop (false,false,nsp,nt,ampl,trace);
	    
	    /* convolve with Ricker wavelet */
	    sf_freqfilt(nt,trace);
	    
	    sf_floatwrite(trace,nt,data);
	}
    }

    exit(0);
}
