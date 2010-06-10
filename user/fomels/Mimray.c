/* 2-D image ray tracing using HWT */
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

int main(int argc, char* argv[])
{
    int nt, nx, nz, it, ix, iy, iz, order;
    float dt, dx, dz, y, z, v, q, x1,x2, z1,z2, r1,r2, xd,zd,rd, d2, a,b, x0;
    float **vv, *xx, *xp, *zz, *zp, *vd, *qq, *rr;
    sf_eno2 vmap;
    sf_file vel, dix;

    sf_init(argc,argv);

    vel = sf_input("in");
    dix = sf_output("out");

    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"o1",&x0)) sf_error("No d1= in input");

    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dz)) sf_error("No d2= in input");
 
    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    
    sf_putint(dix,"n2",4);

    sf_putint(dix,"n3",nt);
    sf_putfloat(dix,"d3",dt);
    sf_putfloat(dix,"o3",0.);
    sf_putstring(dix,"label3","Time");
    sf_putstring(dix,"unit3","s");

    dt *= 0.5; /* one-way time */

    vv = sf_floatalloc2(nx,nz);
    sf_floatread(vv[0],nz*nx,vel);

    if (!sf_getint("order",&order)) order=3;
    /* interpolation order */

    vmap = sf_eno2_init(order,nx,nz);
    sf_eno2_set(vmap,vv);

    xx = sf_floatalloc(nx);
    xp = sf_floatalloc(nx);

    zz = sf_floatalloc(nx);
    zp = sf_floatalloc(nx);
  
    vd = sf_floatalloc(nx);
    qq = sf_floatalloc(nx);
    rr = sf_floatalloc(nx);

    /* on the surface */
    for (ix=0; ix < nx; ix++) {
	zp[ix] = 0.;
	xp[ix] = x0+ix*dx;
	vd[ix] = vv[0][ix];
	qq[ix] = 1.0;
	rr[ix] = vv[0][ix]*dt;
    }

    sf_floatwrite(zp,nx,dix);
    sf_floatwrite(xp,nx,dix);
    sf_floatwrite(vd,nx,dix);
    sf_floatwrite(qq,nx,dix);

    for (it=1; it < nt; it++) {
	/* HWT update */
	for (ix=0; ix < nx; ix++) {
	    if (ix==0) {
		z1 = zp[ix];
		x1 = xp[ix];
		r1 = rr[ix];
	    } else {
		z1 = zp[ix-1];
		x1 = xp[ix-1];
		r1 = rr[ix-1];
	    }
	    if (ix==nx-1) {
		z2 = zp[ix];
		x2 = xp[ix];
		r2 = rr[ix];
	    } else {
		z2 = zp[ix+1];
		x2 = xp[ix+1];
		r2 = rr[ix+1];
	    }
	    
	    zd = z2-z1;
	    xd = x2-x1;
	    rd = r2-r1;
	    d2 = xd*xd+zd*zd;

	    a = rd/d2;
	    b = SF_SIG(xd)*sqrtf(d2-rd*rd)/d2;

	    z1 = a*zd-b*xd;
	    z2 = a*zd+b*xd;

	    x1 = a*xd+b*zd;
	    x2 = a*xd-b*zd;

	    if (z2 < z1) {
		zz[ix] = zp[ix]-rr[ix]*z2;
		xx[ix] = xp[ix]-rr[ix]*x2;
	    } else {
		zz[ix] = zp[ix]-rr[ix]*z1;
		xx[ix] = xp[ix]-rr[ix]*x1;
	    }
	}
	
	for (ix=0; ix < nx; ix++) {
	    zp[ix] = zz[ix];
	    xp[ix] = xx[ix];

	    /* get velocity */
	    y = (xx[ix]-x0)/dx; iy = floorf(y); y -= iy;
	    z = zz[ix]/dz;      iz = floorf(z); z -= iz;

	    if (iy <= 0) {
		iy=0; y=0.0;
	    } else if (iy >= nx) {
		iy=nx-1; y=0.0;
	    } else if (iy == nx-1) {
		y=0.0;
	    }

	    if (iz <= 0) {
		iz=0; z=0.0;
	    } else if (iz >= nz) {
		iz=nz-1; z=0.0;
	    } else if (iz == nz-1) {
		z=0.0;
	    }

	    sf_eno2_apply (vmap,iy,iz,y,z,&v,NULL,FUNC);

	    rr[ix] = v*dt;

	    /* get geometrical spreading */
	    if (ix==0) {
		q = hypotf(xx[ix+1]-xx[ix],
			   zz[ix+1]-zz[ix]);
	    } else if (ix==nx-1) {
		q = hypotf(xx[ix]-xx[ix-1],
			   zz[ix]-zz[ix-1]);
	    } else {
		q = 0.5*hypotf(xx[ix+1]-xx[ix-1],
			       zz[ix+1]-zz[ix-1]);
	    }

	    /* dix velocity */
	    vd[ix] = v;
	    qq[ix] = q/dx;
	}

	sf_floatwrite(zp,nx,dix);
	sf_floatwrite(xp,nx,dix);
	sf_floatwrite(vd,nx,dix);
	sf_floatwrite(qq,nx,dix);
    }
    
    exit(0);
}
