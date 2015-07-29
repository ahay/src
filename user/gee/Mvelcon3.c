/* 3-D finite-difference velocity continuation on a helix */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

#include "velcon3.h"

int main(int argc, char* argv[])
{
    bool adj;
    int inv, nx, ny, nt, nv, i, n;
    float t0, dt, dx, dy, vel;
    float **dat, **img;
    sf_file inp, out;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&ny)) sf_error("No n3= in input");
    n = nt*nx*ny;

    if (!sf_histfloat(inp,"o1",&t0)) t0=0.;
    if (!sf_histfloat(inp,"d1",&dt)) dt=1.;
    if (!sf_histfloat(inp,"d2",&dx)) dx=1.;
    if (!sf_histfloat(inp,"d3",&dy)) dy=dx;

    if (dy != dx) sf_error("dy != dx: %g != %g",dy,dx);

    if (!sf_getbool("adj",&adj)) adj=true;
    /* forward or backward continuation */

    if (!sf_getint("inv",&inv)) inv=1;
    /* inversion type */

    if (!sf_getint("nv",&nv)) nv=nt;
    /* velocity steps */

    if (!sf_getfloat("vel",&vel)) vel=1.;
    /* initial velocity */

    dat = sf_floatalloc2(nt,nx*ny);
    img = sf_floatalloc2(nt,nx*ny);

    if (adj) {
	sf_floatread (dat[0],n,inp); 
	for (i=0; i < n; i++) {
	    img[0][i] = 0.;
	}
    } else {
	sf_floatread (img[0],n,inp); 
	for (i=0; i < n; i++) {
	    dat[0][i] = 0.;
	}
    }

    velcon3_init (inv,vel,0.,t0,nt,nx,ny,nv, dt,dx);
    velcon3_apply (adj,img,dat);

    if (adj) {
	sf_floatwrite (img[0],n,out);
    } else {
	sf_floatwrite (dat[0],n,out);
    }

    exit(0);
}
