/* Constant-velocity prestack exploding reflector in offset. */
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
    int it, nt, ix, nx, iz, nz, ih, nh, ih1, ih2, ihs, snap;
    float dt, dx, dz, dh, v, kx, kz, w, x, c;
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

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1=");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1=");

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2=");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2=");

	if (!sf_histint(data,"n3",&nh)) sf_error("No n3=");
	if (!sf_histfloat(data,"d3",&dh)) sf_error("No d3=");

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* depth samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n3",1); /* stack for now */
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1=");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1=");

	if (!sf_histint(image,"n2",&nx)) sf_error("No n2=");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2=");

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
        /* offset samples (if modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	/* offset sampling (if modeling) */

	sf_putint(data,"n1",nt);
	sf_putfloat(data,"d1",dt);
	sf_putstring(data,"label1","Time");
	sf_putstring(data,"unit1","s");

	sf_putint(data,"n3",nh);
	sf_putfloat(data,"d3",dh);
	sf_putstring(data,"label3","Half-Offset");
	
    }

    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nt,nx);

    prev = sf_floatalloc3(nt,nx,nz);
    curr = sf_floatalloc3(nt,nx,nz);

    if (!sf_getfloat("v",&v)) sf_error("Need v=");
    /* velocity */

    dx = cosft_dk(nx,dx);
    dz = cosft_dk(nz,dz);
    dt = cosft_dk(nt,dt);

   if (NULL != snaps) {
	sf_putint(snaps,"n1",nt);
	sf_putfloat(snaps,"d1",dt);
	sf_putstring(snaps,"label1","Time");

	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putstring(snaps,"label2","Midpoint");

	sf_putint(snaps,"n3",nz);
	sf_putfloat(snaps,"d3",dz);
	sf_putstring(snaps,"label3","Depth");

	sf_putint(snaps,"n4",nh/snap);
	sf_putfloat(snaps,"d4",dh*snap);
	sf_putfloat(snaps,"o4",0.);
	sf_putstring(snaps,"label4","Half-Offset");
    }

   dh *= 2*SF_PI;

    for (iz=0; iz < nz; iz++) {
	for (ix=0; ix < nx; ix++) {
	    for (it=0; it < nt; it++) {
		prev[iz][ix][it] = 0.;
		curr[iz][ix][it] = 0.;
	    }
	}
    }

    if (mig) { /* migration */
	/* initialize image */
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = 0.;
	    }
	}

	/* step backward in offset */
	ih1 = nh-1;
	ih2 = -1;
	ihs = -1;

    } else { /* modeling */
	sf_floatread(img[0],nz*nx,image);

	/* Initialize model */
	for (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		curr[iz][ix][0] = img[ix][iz];
	    }
	}
	cosft3(false,nt,nx,nz,curr);
	
	/* step forward in offset */
	ih1 = 0;
	ih2 = nh;
	ihs = +1;
    }


    /* offset stepping */
    for (ih=ih1; ih != ih2; ih += ihs) {
	sf_warning("ih=%d;",ih);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nt*nx,data);
	    cosft2(false,nt,nx,dat);	    
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (it=0; it < nt; it++) {
		    dat[ix][it] = 0.;
		}
	    }
	}

	if (NULL != snaps && 0 == ih%snap) 
	    sf_floatwrite(curr[0][0],nt*nx*nz,snaps);

	for (iz=0; iz < nz; iz++) {
	    kz = iz*dz;
	    for (ix=0; ix < nx; ix++) {
		kx = (iz==0 && ix==0)? dx: ix*dx;
		x = 4.0f/((kz*kz+kx*kx)*v*v);
		for (it=0; it < nt; it++) {
		    w = it*dt;
		    w = w*w*x-1.0f;

		    c = curr[iz][ix][it];

		    if (mig) {
			c += dat[ix][it];
		    } else {
			dat[ix][it] += c;
		    }

		    if (w >= 0.0f) {
			curr[iz][ix][it] = 2*cosf(kz*dh*sqrtf(w))*c - prev[iz][ix][it];
		    } else {
			curr[iz][ix][it] = c;
		    }
		    prev[iz][ix][it] = c;
		}
	    }
	}

	if (!mig) { /* modeling -> write out data */
	    cosft2(true,nt,nx,dat);
	    sf_floatwrite(dat[0],nt*nx,data);
	}
    }
    sf_warning(".");

    if (mig) {
	for (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (it=0; it < nt; it++) {
		    c = curr[iz][ix][it];
		    img[ix][iz] += c;
		}
	    }
	}

	cosft2(true,nz,nx,img);
	sf_floatwrite(img[0],nz*nx,image);
    }

    exit(0);
}
