/* Natural neighbor interpolation (2-D) */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include<rsf.h>

#include "distance.h"
#include "ntrianglen.h"

int main (int argc,char* argv[]) 
{
    int i, ip, n1, n2, np, ndim, order, n123, *pp, n[2], box[2], *shift[2], b;
    float o1, o2, d1, d2, slow, *dd, **pts, *vv, *h, *bin, *vor, d[2], *rect[2];
    bool isvel;
    sf_upgrad upg;
    sf_file coord, ord, grid, vel;

    sf_init (argc, argv);
    ord = sf_input("in");
    coord = sf_input("coord");
    grid = sf_output("out");

    if (NULL != sf_getstring("velocity")) {
	vel = sf_input("velocity");

	if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in vel");
	if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in vel");
	/* dimensions */
	
	if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in vel");
	if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in vel");
	/* sampling */
	
	if(!sf_histfloat(vel,"o1",&o1)) o1=0.;
	if(!sf_histfloat(vel,"o2",&o2)) o2=0.;
	/* origin */
    } else {
	vel = NULL;

	if(!sf_getint("n1",&n1)) sf_error("Need n1=");
	if(!sf_getint("n2",&n2)) sf_error("Need n2=");
	/* dimensions */
	
	if(!sf_getfloat("d1",&d1)) sf_error("Need d1=");
	if(!sf_getfloat("d2",&d2)) sf_error("Need d2=");
	/* sampling */
	
	if(!sf_getfloat("o1",&o1)) o1=0.;
	if(!sf_getfloat("o2",&o2)) o2=0.;
	/* origin */
    }

    sf_putint(grid,"n1",n1);
    sf_putint(grid,"n2",n2);
    sf_putfloat(grid,"d1",d1);
    sf_putfloat(grid,"d2",d2);
    sf_putfloat(grid,"o1",o1);
    sf_putfloat(grid,"o2",o2);

    n[0]=n1; d[0]=d1;
    n[1]=n2; d[1]=d2;

    if(!sf_getint("order",&order)) order=2;
    /* [1,2] Accuracy order for distance calculation */

    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* if y, the input is velocity; n, slowness squared */

    if (SF_FLOAT != sf_gettype(coord)) sf_error("Need float input");
    if(!sf_histint(coord,"n2",&np)) sf_error("No n2= in input");
    if(!sf_histint(coord,"n1",&ndim) || ndim > 3)  sf_error("Need n1 <= 3 in input");

    pts = sf_floatalloc2 (3,np);
    for (ip=0; ip < np; ip++) {
	sf_floatread(pts[ip],ndim,coord);
	pts[ip][2] = 0.0f;
    }
    
    n123 = n1*n2;

    dd = sf_floatalloc (n123);
    vv = sf_floatalloc (n123);
    pp = sf_intalloc (n123);

    if (NULL != vel) {
	sf_floatread(vv,n123,vel);
	sf_fileclose(vel);

	/* transform velocity to slowness squared */
	if (isvel) {
	    for(i = 0; i < n123; i++) {
		slow = vv[i];
		vv[i] = 1./(slow*slow);
	    }
	} 
    } else {
	for(i = 0; i < n123; i++) {
	    vv[i] = 1.;
	}
    }
    
    /* 1. find distance */
    distance_init (1,n2,n1);  
    distance(np,pts,dd,vv,pp,
	     1,n2,n1,
	     0.,o2,o1,
	     1.,d2,d1,
	     order);

    /* 2. binning */
    sf_int2_init (pts, o1,o2,d1,d2,n1,n2, sf_lin_int, 2, np);
    h = sf_floatalloc(np);

    for (ip=0; ip<np; ip++) {
	h[ip]=1.0f;
    }
    sf_int2_lop (true,false,n123,np,vv,h);

    if (SF_FLOAT != sf_gettype(ord)) sf_error("Need float input");
    sf_floatread(h,np,ord);

    bin = sf_floatalloc(n123);
    
    sf_int2_lop (true,false,n123,np,bin,h);
    
    box[0] = 1;
    for (i=0; i < n123; i++) {
	/* normalize by the fold */
	if (vv[i] > FLT_EPSILON) bin[i] /=vv[i];
	vv[i]=0.0f;
	pp[i]=0;
	b = ceilf(dd[i]);
	if (b > box[0]) box[0] = b;
    }
    box[1] = box[0];

    /* 3. voronoi interpolation */

    vor = sf_floatalloc(n123);

    upg = sf_upgrad_init(2,n,d);

    sf_upgrad_set(upg,dd);
    sf_upgrad_solve(upg,vv,vor,bin); 
	
    /* 4. smoothing */

    rect[0] = dd; shift[0] = pp;
    rect[1] = dd; shift[1] = pp;

    ntrianglen_init(2,box,n,rect,shift,1);
    ntrianglen_lop(false,false,n123,n123,vor,bin);
	    
    sf_floatwrite(bin,n123,grid); 

    exit (0);
}

