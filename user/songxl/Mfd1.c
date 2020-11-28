/* 1-D Optimized finite-difference wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
    int nx, nt, ix, it, isx;
    float dt, dx;
    float *old, *nxt, *cur, *sig, *v;
    sf_file in, out, vel;
    int im,im2,ip,ip2;

    sf_init(argc,argv);
    in  = sf_input("in");
    vel = sf_input("vel");   /* velocity */
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_getint("isx",&isx)) isx=(int)(nx/2);
    if (!sf_getint("nt",&nt)) sf_error("No nt in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt in input");

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 

    sig = sf_floatalloc(nx);
    old = sf_floatalloc(nx);
    nxt = sf_floatalloc(nx);
    cur = sf_floatalloc(nx);
    v   = sf_floatalloc(nx);
    
 
    sf_floatread(v,nx,vel);
    sf_floatread(sig,nx,in);		

	/* initial conditions */
    for (ix=0; ix < nx; ix++){
        cur[ix] =  sig[ix];
        old[ix] =  0.0; 
	nxt[ix] = 0.;
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
	sf_floatwrite(cur,nx,out);
	/* Stencil */
	for (ix=0; ix < nx; ix++) {
            im  =(ix-1+nx)%nx; 
            im2 =(ix-2+nx)%nx; 
            ip = (ix+1+nx)%nx;
            ip2 =(ix+2+nx)%nx;
	    nxt[ix] = dt*dt/(12.0*dx*dx)*(-30.0*cur[ix] +16.0* (cur[im]+cur[ip]) - (cur[im2]+cur[ip2]))*v[ix]*v[ix] 
                      +2.0*cur[ix]- old[ix];
	}
       for (ix=0; ix < nx; ix++) {
	    old[ix] = cur[ix];
	    cur[ix] = nxt[ix];
	}
    }
    exit(0);
}
