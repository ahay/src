/* Computing distance function by fast marching eikonal solver (3-D). */
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

#include<rsf.h>

#include "distance.h"

int main (int argc,char* argv[]) 
{
    int i, n1, n2, n3, np, ndim, order, n123, *pp;
    float o1, o2, o3, d1, d2, d3, slow, *dd, **pts, *vv;
    bool isvel;
    sf_file points, dist, vel;

    sf_init (argc, argv);
    points = sf_input("in");
    dist = sf_output("out");

    if (NULL != sf_getstring("velocity")) {
	vel = sf_input("velocity");

	if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in vel");
	if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in vel");
	if(!sf_histint(vel,"n3",&n3)) n3=1;
	/* dimensions */
	
	if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in vel");
	if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in vel");
	if(!sf_histfloat(vel,"d3",&d3)) d3=d2;
	/* sampling */
	
	if(!sf_histfloat(vel,"o1",&o1)) o1=0.;
	if(!sf_histfloat(vel,"o2",&o2)) o2=0.;
	if(!sf_histfloat(vel,"o3",&o3)) o3=0.;
	/* origin */
    } else {
	vel = NULL;

	if(!sf_getint("n1",&n1)) sf_error("Need n1=");
	if(!sf_getint("n2",&n2)) sf_error("Need n2=");
	if(!sf_getint("n3",&n3)) n3=1;
	/* dimensions */
	
	if(!sf_getfloat("d1",&d1)) sf_error("Need d1=");
	if(!sf_getfloat("d2",&d2)) sf_error("Need d2=");
	if(!sf_getfloat("d3",&d3)) d3=d2;
	/* sampling */
	
	if(!sf_getfloat("o1",&o1)) o1=0.;
	if(!sf_getfloat("o2",&o2)) o2=0.;
	if(!sf_getfloat("o3",&o3)) o3=0.;
	/* origin */
    }

    sf_putint(dist,"n1",n1);
    sf_putint(dist,"n2",n2);
    sf_putint(dist,"n3",n3);
    sf_putfloat(dist,"d1",d1);
    sf_putfloat(dist,"d2",d2);
    sf_putfloat(dist,"d3",d3);
    sf_putfloat(dist,"o1",o1);
    sf_putfloat(dist,"o2",o2);
    sf_putfloat(dist,"o3",o3);

    if(!sf_getint("order",&order)) order=2;
    /* [1,2] Accuracy order */

    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* if y, the input is velocity; n, slowness squared */

    if (SF_FLOAT != sf_gettype(points)) sf_error("Need float input");
    if(!sf_histint(points,"n2",&np)) sf_error("No n2= in input");
    if(!sf_histint(points,"n1",&ndim) || ndim != 3) 
	sf_error("Need n1=3 in input");

    pts = sf_floatalloc2 (ndim,np);
    sf_floatread(pts[0],np*ndim,points);
    
    n123 = n1*n2*n3;

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
    
    distance_init (n3,n2,n1,np);
  
    distance(np,pts,dd,vv,pp,
	     n3,n2,n1,
	     o3,o2,o1,
	     d3,d2,d1,
	     order);
	
    sf_floatwrite (dd,n123,dist);


    exit (0);
}

/* 	$Id: Meikonal.c 825 2004-10-07 08:11:17Z fomels $	 */
