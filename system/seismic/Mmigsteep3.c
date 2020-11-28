/* 3-D Kirchhoff time migration for antialiased steep dips. 

Combine with sfmig3 antialias=flat for the complete response.
*/
/*
  Copyright (C) 2007 University of Texas at Austin

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
    int nt, nx, ny, ix, iy, n1, it, iz, ntr, itr;
    float hdr[3], *trace=NULL, *x1=NULL, *x2=NULL, *amp=NULL, *xmod=NULL, *ymod=NULL, ***out=NULL;
    float h_in, x_in, y_in, vel, vel2;
    float x,y, z, t, dx, dy, ox, oy, ry, h, sq, tx, ot, dt, ty, rx;
    sf_map xstr, ystr;
    sf_file in=NULL, mig=NULL, head=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    head = sf_input("hdr");
    mig = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    ntr = sf_leftsize(in,1);

    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");

    if (!sf_getint("n2",&nx))   sf_error("Need n2="); 
    if (!sf_getfloat("d2",&dx)) sf_error("Need d2="); 
    if (!sf_getfloat("o2",&ox)) sf_error("Need o2="); 
	sf_putint(mig,"n2",nx);
	sf_putfloat(mig,"d2",dx);
	sf_putfloat(mig,"o2",ox);

    if (!sf_getint("n3",&ny))   sf_error("Need n3="); 
    if (!sf_getfloat("d3",&dy)) sf_error("Need d3="); 
    if (!sf_getfloat("o3",&oy)) sf_error("Need o3="); 
	sf_putint(mig,"n3",ny);
	sf_putfloat(mig,"d3",dy);
	sf_putfloat(mig,"o3",oy);

    if (!sf_getint("n1",&n1))   sf_error("Need n1="); 
	sf_putint(mig,"n1",n1);

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* migration velocity */
    vel = 2./vel; 
    vel2 = vel*vel;

    trace = sf_floatalloc(nt);
    x1 = sf_floatalloc(nt);
    x2 = sf_floatalloc(nt);
    amp = sf_floatalloc(nt);
    xmod = sf_floatalloc(nx);
    ymod = sf_floatalloc(ny);
    out = sf_floatalloc3(n1,nx,ny);

    for (it=0; it < n1*nx*ny; it++) {
	out[0][0][it] = 0.;
    }

    xstr = sf_stretch_init (nx,ox,dx,nt,0.01,true);
    ystr = sf_stretch_init (ny,oy,dy,nt,0.01,true);

    for (itr=0; itr < ntr; itr++) {
	sf_floatread (trace,nt,in);
	sf_floatread (hdr,3,head);

	x_in = hdr[0];
	y_in = hdr[1];
	h_in = hdr[2];
	
	h = h_in * h_in * vel2;

	for (iz=0; iz < n1; iz++) {
	    z = ot + iz*dt;
	    z *= z;

	    for (iy=0; iy < ny; iy++) {
		y = oy + iy*dy - y_in;
		ry = fabsf(y*dy) * vel2;
		y = y*y * vel2 + z;
		
		for (it=0; it < nt; it++) {
		    x1[it] = ox - 10.*dx;
		    x2[it] = ox - 10.*dx;
		    amp[it] = 0.;
		}
		for (it=0; it < nt; it++) {
		    t = ot + it*dt;
		    sq = t*t - h; 
		    if (sq <= y) continue;

		    x = t*sqrtf((sq - y)/sq);
		    ty = t*dt*(sq*sq+ h*y)/(sq*sq*x);
		    tx = SF_MAX(ty, t*ry/sqrtf(sq*(sq - y)));

		    if (tx >= dx*vel) continue;

		    x1[it] = x_in+x/vel;
		    x2[it] = x_in-x/vel;
		    amp[it] = trace[it]*dx*vel/ty;
		}
		sf_stretch_define (xstr,x1);
		sf_stretch_apply (xstr,amp,xmod);
		for (ix=0; ix < nx; ix++) {
		    out[iy][ix][iz] += xmod[ix];
		}
		sf_stretch_define (xstr,x2);
		sf_stretch_apply (xstr,amp,xmod);
		for (ix=0; ix < nx; ix++) {
		    out[iy][ix][iz] += xmod[ix];
		}
	    } /* iy */
	    for (ix=0; ix < nx; ix++) {
		x = ox + ix*dx - x_in;
		rx = fabsf(x*dx) * vel2;
		x *= x * vel2;

		for (it=0; it < nt; it++) {
		    x1[it] = oy - 10.*dy;
		    x2[it] = oy - 10.*dy;
		    amp[it] = 0.;
		}

		for (it=0; it < nt; it++) {
		    t = ot + it*dt;
		    ry = fabsf(t*dt);
		    t *= t; 
		    if (t <= 0.) continue;
		    sq = t + x*h/t - h - z - x; 
		    if (sq <= 0.) continue;

		    y = sqrtf(sq);
		    ty = ry*fabsf(1.0f-x*h/(t*t))/y;
		    tx = SF_MAX(ty, rx*fabsf(1.0f-h/t)/y);

		    if (tx >= dy*vel) continue;

		    x1[it] = y_in+y/vel;
		    x2[it] = y_in-y/vel;
		    amp[it] = trace[it]*dy*vel/ty;
		}
		sf_stretch_define (ystr,x1);
		sf_stretch_apply (ystr,amp,ymod);
		for (iy=0; iy < ny; iy++) {
		    out[iy][ix][iz] += ymod[iy];
		}
		sf_stretch_define (ystr,x2);
		sf_stretch_apply (ystr,amp,ymod);
		for (iy=0; iy < ny; iy++) {
		    out[iy][ix][iz] += ymod[iy];
		}
	    } /* ix */
	} /* iz */
    } /* itr */

    sf_floatwrite(out[0][0],n1*nx*ny,mig);

    exit(0);
}

