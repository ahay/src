/* Convert Dix velocity to interval velocity (2-D). */
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

#include "cameron.h"

int main(int argc, char *argv[])
{
    int ix, it, i, i1, nx, nt, nxt, niter;
    bool *k;
    float x0, dx, dt, *v, *x, *z, *r;
    sf_file vdix, vint, xz;

    sf_init(argc,argv);
    vdix = sf_input("in");
    vint = sf_output("out");

    xz = sf_output("xz");
    sf_putint(xz,"n3",2);

    if (!sf_histint(vdix,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(vdix,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(vdix,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(vdix,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vdix,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("niter",&niter)) niter=nt*3;
    /* number of iterations */

    nxt = nx*nt;

    v = sf_floatalloc(nxt);
    x = sf_floatalloc(nxt);
    z = sf_floatalloc(nxt);
    r = sf_floatalloc(2*nxt);
    k = sf_boolalloc(nxt);

    sf_floatread(v,nxt,vdix);

    cameron_init(nt,nx,dt,dx,v);

    for (ix=0; ix < nx; ix++) {
	i = ix*nt;

	x[i]   = x0+ix*dx;
	x[i+1] = x0+ix*dx;
	z[i]   = 0.;
	z[i+1] = dt*v[i];
	k[i]   = true;
	k[i+1] = true;

	for (it=2; it < nt; it++) {
	    i1 = i+it;

	    x[i1] = x[i];
	    z[i1] = z[i1-1]+dt*v[i1-1];
	    k[i1] = (ix==0) || (ix==nx-1);
	}
    }

    for (i=0; i < 2*nxt; i++) {			
	r[i] = 0.;
    }

    sf_solver (cameron_lop,sf_cgstep,nxt,2*nxt,x,r,niter,"x0",x,"known",k,"verb",true,"end");
    sf_cgstep_close();
    sf_solver (cameron_lop,sf_cgstep,nxt,2*nxt,z,r,niter,"x0",z,"known",k,"verb",true,"end");
    sf_cgstep_close();

    sf_floatwrite(x,nxt,xz);
    sf_floatwrite(z,nxt,xz);

    for (ix=0; ix < nx; ix++) {
	i = ix*nt;
	for (it=1; it < nt; it++) {
	    i1=i+it;

	    v[i1] = hypotf(x[i1]-x[i1-1],
			   z[i1]-z[i1-1])/dt;
	}
    }

    sf_floatwrite(v,nxt,vint);
    exit(0);
}
