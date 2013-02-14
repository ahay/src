/* 2-D multiple arrivals by cell ray tracing. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "agrid.h"

static float da, a0;
static sf_celltrace ct;

static void raytrace (float q, void* x, float* xzt);

int main(int argc, char* argv[])
{
    int nx, na, na2, ia, nz, order, maxsplit, ix, iz, *siz;
    float **place, *slow, **out, dx,dz, x0,z0, x[2];
    float max1, min1, max2, min2;
    bool isvel, lsint;
    agrid grd;
    sf_file vel, outp, size, grid;
        
    sf_init (argc,argv);

    /* get 2-D grid parameters */
    vel = sf_input("in");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&x0)) x0=0.;

    outp = sf_output("out");
    
    sf_putint(outp,"n4",nz);
    sf_putfloat(outp,"d4",dz);
    sf_putfloat(outp,"o4",z0);

    sf_putint(outp,"n3",nx);
    sf_putfloat(outp,"d3",dx);
    sf_putfloat(outp,"o3",x0);

    if (!sf_getint("na",&na)) na=60;
    /* number of angles */
    if (!sf_getfloat("da",&da)) da=3.1;
    /* angle increment (in degrees) */
    if (!sf_getfloat("a0",&a0)) a0=-90.;
    /* initial angle (in degrees) */

    if (!sf_getint("maxsplit",&maxsplit)) maxsplit=10;
    /* maximum splitting for adaptive grid */

    if (!sf_getfloat("minx",&min1)) min1=0.5*dx;
    /* parameters for adaptive grid */
    if (!sf_getfloat("maxx",&max1)) max1=2.*dx;
    if (!sf_getfloat("mina",&min2)) min2=0.5*da;
    if (!sf_getfloat("maxa",&max2)) max2=2.*da;

    sf_putint(outp,"n2",na);
    sf_putfloat(outp,"d2",da);
    sf_putfloat(outp,"o2",a0);

    da *= (SF_PI/180.);
    a0 *= (SF_PI/180.);

    sf_putint(outp,"n1",5);

    size = sf_output("size");
    sf_putint(size,"n1",nx);
    sf_putint(size,"n2",nz);
    sf_settype(size,SF_INT);

    grid = sf_output("grid");

    /* additional parameters */
    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* y: velocity, n: slowness */
    if(!sf_getint("order",&order)) order=3;
    /* velocity interpolation order */
    if (!sf_getbool("lsint",&lsint)) lsint=false;
    /* if use least-squares interpolation */

    slow  = sf_floatalloc(nz*nx);
    place = sf_floatalloc2(5,na);
    siz  = sf_intalloc(nx);

    sf_floatread(slow,nx*nz,vel);

    if (isvel) {
      for(ix = 0; ix < nx*nz; ix++){
	slow[ix] = 1./slow[ix];
      }
    }

    ct = sf_celltrace_init (lsint, order, nz*nx, nz, nx, dz, dx, z0, x0, slow);
    free (slow);

    grd = agrid_init (na, 5, maxsplit);
    agrid_set (grd,place);

    for (iz = 0; iz < nz; iz++) {
	x[0] = z0+iz*dz;

	sf_warning("depth %d of %d;",iz+1, nz);
	for (ix = 0; ix < nx; ix++) {
	    x[1] = x0+ix*dx;

	    fill_grid(grd,min1,max1,min2,max2,(void*) x, raytrace);
	
	    na2 = grid_size (grd);
	    out = write_grid (grd);

	    siz[ix] = na2;
	
	    for (ia=0; ia < na2; ia++) {
	      if (ia < na)
		  sf_floatwrite (place[ia], 5, outp);

	      sf_floatwrite(out[ia], 5, grid);
	    }
	    free (out);
	}
	
	sf_intwrite (siz,nx,size);
    }
    sf_warning(".");

    exit (0);
}

static void raytrace (float q, void* xv, float* xzt) {
  float x[2], p[2], a, *xf;
    int it;

    xf = (float*) xv;
    x[0] = xf[0];
    x[1] = xf[1];

    a = a0+q*da;
    p[1] = sin(a);
    p[0] = -cos(a);

    xzt[3] = p[1];
    xzt[2] = sf_cell_trace (ct,(float*) x,p,&it,NULL);
    xzt[1] = x[0];
    xzt[0] = x[1];
    xzt[4] = sf_cell_p2a(p)*180./SF_PI;
}

/* 	$Id: Mtrace2.c 853 2004-11-05 12:54:22Z fomels $	 */
