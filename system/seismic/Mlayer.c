/* Ray tracing in a layered medium. */
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
#include <rsf.h>

#include "layer.h"

int main(int argc, char* argv[]) 
{
    int nx, nc, nxc, ix, ic, nv, iv;
    float **rfl, **crv, **dip, *v0, **g, **r0;
    float slow, dx, x0, s0, s1;
    char **type;
    bool twod;
    sf_file curv;
    
    sf_init(argc,argv);
    curv = sf_input("in");
    
    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&nc) || nc < 2) 
	sf_error("Need n2 >= 2 in input");
    /* number of interfaces */
    nxc = nx*nc;

    /*** Initialize shots ***/
    if (!sf_getfloat("s0",&s0)) sf_error("Need s0=");
    /* first point */
    if (!sf_getfloat("s1",&s1)) sf_error("Need s1=");
    /* last point */
    
    /*** Initialize reflector ***/

    rfl = sf_floatalloc2(nx,nc);
    dip = sf_floatalloc2(nx,nc);
    crv = sf_floatalloc2(nx,nc);
    
    sf_floatread(crv[0],nxc,curv);
    sf_fileclose(curv);

    /* reflector dip */
    if (NULL != sf_getstring("dip")) {
	curv = sf_input("dip");
	sf_floatread(dip[0],nxc,curv);
	sf_fileclose(curv);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		dip[ic][ix] = 0.;
	    }
	}
    }

    /* reflector curvature */
    if (NULL != sf_getstring("crv")) {
	curv = sf_input("crv");
	sf_floatread(crv[0],nxc,curv);
	sf_fileclose(curv);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		crv[ic][ix] = 0.;
	    }
	}
    }  
    
    /*** Initialize velocity ***/
    nv = nc-1;
    v0 = sf_floatalloc(nv);
    g = sf_floatalloc2(nv,2);
    r0 = sf_floatalloc2(nv,2);
    type = (char**) sf_alloc(nv,sizeof(char*));
 
    if (!sf_getfloats("vel",v0,nv)) sf_error("Need vel=");
    /* velocity */
    if (!sf_getfloat("gradx",g[0])) { /* velocity gradient in x */
	for (iv=0; iv < nv; iv++) {
	    g[0][iv] = 0.;
	}
    }
    if (!sf_getfloat("gradz",g[1])) { /* velocity gradient in z */
	for (iv=0; iv < nv; iv++) {
	    g[1][iv] = 0.;
	}
    }

    if (!sf_getstrings("type",type,nv)) {
	/* type of velocity ('c': constant [default], 
	   's': linear sloth, 'v': linear velocity) */
	for (iv=0; iv < nv; iv++) {
	    type[iv] = "const";
	}
    } 
    for (iv=0; iv < nv; iv++) {
	if ('s'==type[iv][0]) {
	    /* linear slowness squared */
	    slow = 1./(v0[iv]*v0[iv]);
	    /* slowness squared */
	    g[0][iv] *= -2.*slow/v0[iv];
	    g[1][iv] *= -2.*slow/v0[iv];
	    v0[iv] = slow;
	} 
    }

    if (!sf_getbool("twod",&twod)) twod=false;
    /* 2-D or 2.5-D */
	
    if (!sf_getfloat("refx",r0[0])) { /* velocity reference in x */
	for (iv=0; iv < nv; iv++) {
	    r0[0][iv] = x0;
	}
    }
    if (!sf_getfloat("refz",r0[1])) { /* velocity reference in z */
	for (iv=0; iv < nv; iv++) {
	    r0[1][iv] = 0.;
	}
    }


    exit(0);
}
    
