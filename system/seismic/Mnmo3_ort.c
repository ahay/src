/* 3-D Normal moveout using orthogonal parametrization

input data has gathers along *4th* axis; 
velocity file contains slowness squared with n2=3 (Wavg,Wcos,Wsin);
offset file contains x,y offset pairs for input data
*/
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
#include "fint1.h"

static float *vx, *vy, *vxy, x, y;

/*"NMO-Ellipse" 3D Correction*/
static float nmo_map(float t, int it) {
    float vxi, vyi, vxyi;
    
    vxi = vx[it];
    vyi = vy[it];
    vxyi = vxy[it];

    vxyi = vxi*(x*x+y*y) + vyi*(x*x-y*y) + 2*vxyi*x*y;
    t = t*t + vxyi;
    if (t > 0.) t=sqrtf(t);
    return t;
}

int main (int argc, char* argv[])
{
    fint1 nmo; /* using cubic spline interpolation */
    mapfunc nmofunc;
    bool half;
    int ix,iy, nt, nx, ny, nh,ih, ntr, itr, nw, i4, n4, n, mute;
    float dt, t0, x0, y0, dx, dy, eps, scale;
    float *trace, *offx, *offy;
    sf_file cmp, nmod, vel, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    nmod = sf_output("out");

    nmofunc =  nmo_map;
    
    /* Read-in input dimensions */
    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histint(cmp,"n4",&n4)) sf_warning("Only one gather detected (4th axis).");

    n4 = sf_leftsize(cmp,3);
    ntr=nx*ny*n4;

    /* Read-in velocity dimensions */
    vel = sf_input("velocity");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float velocity");
    if (!sf_histint(vel,"n1",&n) || nt != n) sf_error("Need n1=%d in velocity",nt);
    if (!sf_histint(vel,"n2",&n) || 3 != n) sf_error("Need n2=3 in velocity");
    if (n4 != sf_leftsize(vel,2)) sf_error("Wrong dimensions in velocity");
    

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second and third axes are half-offset instead of full offset */

    /* Read-in or create offset vectors */
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nh = sf_filesize(offset);
	if (nh != 2*ntr) sf_error("Wrong dimensions in offset");
	if (!sf_histint(offset,"n4",&n) || 2 != n) sf_error("Need n4=2 in offset");
	offx = sf_floatalloc(ntr);
	offy = sf_floatalloc(ntr);
	sf_floatread (offx,ntr,offset);
	sf_floatread (offy,ntr,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d3",&dy)) sf_error("No d3= in input");
	if (!sf_histfloat(cmp,"o3",&y0)) sf_error("No o3= in input");

	nh = ntr;
	offx = sf_floatalloc(nh);
	offy = sf_floatalloc(nh);
	ih=0;
	for (i4=0; i4 < n4; i4++) {
	    for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
		    offy[ih] = y0 + iy*dy;
		    offx[ih] = x0 + ix*dx;
		    ih+=1;
		}
	    }
	}
	offset = NULL;
    }

    if (half) { /* get full offset - multiply by 2 */
	scale=2.0;
	sf_warning("Since half=y, offsets were doubled.");
    }else{
	scale=1.0;
      sf_warning("Since half=n, offsets not doubled.");
    }
	
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    vx = sf_floatalloc(nt);
    vy = sf_floatalloc(nt);
    vxy = sf_floatalloc(nt);


    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = fint1_init (nw, nt, mute);    
    ntr=nx*ny;
    ih=0;

    for (i4=0; i4 < n4; i4++) { /* loop over CMPs */

	sf_warning("CMP %d of %d;",i4+1,n4);

	/*Read velocity model for CMP */
	sf_floatread (vx,nt,vel);
	sf_floatread (vy,nt,vel);
	sf_floatread (vxy,nt,vel);
	
	for (itr = 0; itr < ntr; itr++) { /* loop over traces */
	    
	    /*Offsets*/
	    x = scale*offx[ih];
	    y = scale*offy[ih];
	    
            /*Trace*/
	    sf_floatread (trace,nt,cmp);
	    
            /*Apply NMO correction*/
	    fint1_set(nmo,trace); 
	    stretch(nmo,nmofunc,nt,dt,t0,nt,dt,t0,trace,eps);

	    /*Output corrected trace*/
	    sf_floatwrite (trace,nt,nmod);
	    
	    ih+=1;
	    
	} /* itr */
    } /* i4 */
    
    exit (0);
}


