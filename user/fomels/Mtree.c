/* Multiple arrivals with a fast algorithm. */
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

#include <math.h>

#include <rsf.h>

#include "tree.h"

int main(int argc, char* argv[])
{
    bool velocity, debug;
    int nz,nx, i, iz, na, ix, order, naxz;
    float **slow, **node;
    float dz,dx,da,x0,z0,a0;
    sf_file vel, out;

    sf_init (argc,argv);
    vel = sf_input ("in");
    out = sf_output("out");

    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
  
    if (!sf_histfloat(vel,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(vel,"o2",&x0)) sf_error("No o2= in input");

    sf_putint(out,"n4",nz);
    sf_putfloat(out,"d4",dz);
    sf_putfloat(out,"o4",z0);

    sf_putint(out,"n3",nx);
    sf_putfloat(out,"d3",dx);
    sf_putfloat(out,"o3",x0);

    if (!sf_getint("na",&na)) na=361;
    /* number of angles */
    if (!sf_getfloat("da",&da)) da=1.;
    /* angle increment (in degrees) */
    if (!sf_getfloat("a0",&a0)) a0=-180.;
    /* first angle (in degrees) */

    sf_putint(out,"n2",na);
    sf_putfloat(out,"d2",da);
    sf_putfloat(out,"o2",a0);

    sf_putint(out,"n1",4);

    da *= (SF_PI/180.);
    a0 *= (SF_PI/180.);

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* y: theinput is velocity; n: slowness */
    if(!sf_getint("order",&order)) order=3; 
    /* accuracy order */

    if (!sf_getbool("debug",&debug)) debug=false;
    /* debugging flag */

    slow  = sf_floatalloc2(nz,nx);
    sf_floatread(slow[0],nz*nx,vel);

    if (velocity) { /* transform velocity to slowness */
	for(ix = 0; ix < nx; ix++){
	    for (iz = 0; iz < nz; iz++) {
		slow[ix][iz] = 1./slow[ix][iz];
	    }
	}
    }
   
    naxz = na*nx*nz;
    node = sf_floatalloc2(4,naxz);
    tree_init (order, nz, nx, na, 
	       dz, dx, da, z0, x0, a0, slow, node);
    
    for (i = 0; i < naxz; i++) {
	node[i][0]=0.;
	node[i][1]=0.;
	node[i][2]=-1.;
	node[i][3]=0.;
    }

    /* build dependency graph */
    
    tree_build(debug);

/*    tree_print(); */

    sf_floatwrite(node[0],4*naxz,out);
    
    exit (0);
}

/* 	$Id: Mtree.c 1575 2005-11-21 14:09:06Z fomels $	 */
