/* Kjartansson-style tomography */
/*
 Copyright (C) 2014 The University of Texas at Austin
 
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

#include "xtomo.h"

int main(int argc, char* argv[])
{
    bool adj;
    int nt, nz, nx, nh, ny, niter, nm, nd;
    float ot, oz, ox, oh, oy;
    float dt, dz, dx, dh, dy, zmax;
    float *modl, *data;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    zmax = ot;

    if (adj) {
	if (!sf_histint(inp,"n2",&nh)) sf_error("No n2= in input");
	if (!sf_histfloat(inp,"o2",&oh)) sf_error("No o2= in input");
	if (!sf_histfloat(inp,"d2",&dh)) sf_error("No d2= in input");

	if (!sf_histint(inp,"n3",&ny)) sf_error("No n3= in input");
	if (!sf_histfloat(inp,"o3",&oy)) sf_error("No o3= in input");
	if (!sf_histfloat(inp,"d3",&dy)) sf_error("No d3= in input");

	if (!sf_getint("nz",&nz))   nz=nh;
	if (!sf_getfloat("oz",&oz)) oz=oh;
	if (!sf_getfloat("dz",&dz)) dz=dh;

	if (!sf_getint("nx",&nx))   nx=ny;
	if (!sf_getfloat("ox",&ox)) ox=oy;
	if (!sf_getfloat("dx",&dx)) dx=dy;
    } else {
	if (!sf_histint(inp,"n2",&nz)) sf_error("No n2= in input");
	if (!sf_histfloat(inp,"o2",&oz)) oz=-zmax;
	if (!sf_histfloat(inp,"d2",&dz)) dz=(zmax-oz)/(nz-1);

	if (!sf_histint(inp,"n3",&nx)) sf_error("No n3= in input");
	if (!sf_histfloat(inp,"o3",&ox)) sf_error("No o3= in input");
	if (!sf_histfloat(inp,"d3",&dx)) sf_error("No d3= in input");
	
	if (!sf_getint("nh",&nh))   nh=nz;
	if (!sf_getfloat("oh",&oh)) oh=oz;
	if (!sf_getfloat("dh",&dh)) dh=dz;

	if (!sf_getint("ny",&ny))   ny=nx;
	if (!sf_getfloat("oy",&oy)) oy=ox;
	if (!sf_getfloat("dy",&dy)) dy=dx;
    }

    nm = nt*nz*nx;
    nd = nt*nh*ny;

    modl = sf_floatalloc(nm);
    data = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=-1;
    /* number of iterations */

    oh /= 2;
    dh /= 2;
    /* external data representation is full offset */

    if( oh < 0.) { /*  Assume off-end data is at pos offset */ 
	oh= -oh;  
	dh= -dh;
    }    
/*    hmax = SF_MAX(oh, oh+dh*(nh-1)); */

    if (adj) {
	sf_floatread(data,nd,inp);
    } else {
	sf_floatread(modl,nm,inp);
    }
    
    xtomo_init(oz,dz, 
	       ox,dx, 
	       oh,dh,
	       oy,dy,
	       nt,nz,nx,nh,ny);

    

    if (niter < 0 ) {
	xtomo_lop(adj,false,nm,nd,modl,data);
    } else {
	sf_solver (xtomo_lop, sf_cgstep, nm, nd, modl, data, 
		   niter, "verb", true, "end"); 
    }

    if (adj) {
	sf_floatwrite(modl,nm,out);
    } else {
	sf_floatwrite(data,nd,out);
    }

    exit(0);
}
