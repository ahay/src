/* Constant-velocity prestack exploditing reflector. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "cosft3.h"

int main(int argc, char* argv[])
{
    bool mig;
    int it, nt, ix, nx, iz, nz, ih, nh, it1, it2, its, snap;
    float dt, dx, dz, dh, v, kx, kz, kh, h, x, c;
    float ***prev, ***curr, **img, **dat;
    sf_file data, image, snaps;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */

    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    
    snaps = (snap > 0)? sf_output("snaps"): NULL;
    /* (optional) snapshot file */

    if (mig) { /* migration */
	data = sf_input("in");
	image = sf_output("out");

	if (!sf_histint(data,"n1",&nh));
	if (!sf_histfloat(data,"d1",&dh));

	if (!sf_histint(data,"n2",&nx));
	if (!sf_histfloat(data,"d2",&dx));

	if (!sf_histint(data,"n3",&nt));
	if (!sf_histfloat(data,"d3",&dt));

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* time samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* time sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n3",1); /* stack for now */
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_histint(image,"n1",&nz));
	if (!sf_histfloat(image,"d1",&dz));

	if (!sf_histint(image,"n2",&nx));
	if (!sf_histfloat(image,"d2",&dx));

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
        /* offset samples (if modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	/* offset sampling (if modeling) */

	sf_putint(data,"n1",nh);
	sf_putfloat(data,"d1",dh);
	sf_putstring(data,"label1","Half-Offset");

	sf_putint(data,"n3",nt);
	sf_putfloat(data,"d3",dt);
	sf_putstring(data,"label3","Time");
	sf_putstring(data,"unit3","s");
    }

    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nh,nx);

    prev = sf_floatalloc3(nh,nx,nz);
    curr = sf_floatalloc3(nh,nx,nz);

    if (!sf_getfloat("v",&v)) sf_error("Need v=");
    /* velocity */

    v *= SF_PI*dt;

    dx = cosft_dk(nx,dx);
    dz = cosft_dk(nz,dz);
    dh = cosft_dk(nh,dh);

   if (NULL != snaps) {
	sf_putint(snaps,"n1",nh);
	sf_putfloat(snaps,"d1",dh);
	sf_putstring(snaps,"label1","Half-Offset");

	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putstring(snaps,"label2","Midpoint");

	sf_putint(snaps,"n3",nz);
	sf_putfloat(snaps,"d3",dz);
	sf_putstring(snaps,"label3","Depth");

	sf_putint(snaps,"n4",nt/snap);
	sf_putfloat(snaps,"d4",dt*snap);
	sf_putfloat(snaps,"o4",0.);
	sf_putstring(snaps,"label4","Time");
    }

    for (iz=0; iz < nz; iz++) {
	for (ix=0; ix < nx; ix++) {
	    for (ih=0; ih < nh; ih++) {
		prev[iz][ix][ih] = 0.;
		curr[iz][ix][ih] = 0.;
	    }
	}
    }

    if (mig) { /* migration */
	/* initialize image */
	for (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		img[ix][iz] = 0.;
	    }
	}

	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;

    } else { /* modeling */
	sf_floatread(img[0],nz*nx,image);

	/* Initialize model */
	for (iz=1; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		curr[iz][ix][0] = img[ix][iz];
	    }
	}
	cosft3(false,nh,nx,nz,curr);
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }


    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nh*nx,data);
	    cosft2(false,nh,nx,dat);	    
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    dat[ix][ih] = 0.;
		}
	    }
	}

	if (NULL != snaps && 0 == it%snap) 
	    sf_floatwrite(curr[0][0],nh*nx*nz,snaps);

	for (iz=1; iz < nz; iz++) {
	    kz = iz*dz;
	    for (ix=0; ix < nx; ix++) {
		kx = ix*dx;
		x = (kz*kz+kx*kx)/kz;
		for (ih=0; ih < nh; ih++) {
		    kh = ih*dh;
		    h = (kz*kz+kh*kh)/kz;
		    
		    c = curr[iz][ix][ih];

		    if (mig) {
			c += (iz==nz-1)? dat[ix][ih]*0.5: dat[ix][ih];
		    } else {
			dat[ix][ih] += (iz==nz-1)? c*0.5: c;
		    }

		    curr[iz][ix][ih] = 2*cosf(v*sqrtf(x*h))*c - prev[iz][ix][ih];
		    prev[iz][ix][ih] = c;
		}
	    }
	}

	if (!mig) { /* modeling -> write out data */
	    cosft2(true,nh,nx,dat);
	    sf_floatwrite(dat[0],nh*nx,data);
	}
    }
    sf_warning(".");

    if (mig) {
	for (iz=1; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    c = curr[iz][ix][ih];
		    img[ix][iz] += (iz==nz-1)? c*0.5: c;
		}
	    }
	}

	cosft2(true,nz,nx,img);
	sf_floatwrite(img[0],nz*nx,image);
    }

    exit(0);
}
