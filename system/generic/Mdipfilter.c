/* Filter data based on dip in 2-D or 3-D.

February 2014 program of the month:
http://ahay.org/blog/2014/02/06/program-of-the-month-sfdipfilter/
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

#include <float.h>
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nw, nx, ny, iw, ix, iy, i3, n3, dim;
    float x0,dx,x, y0,dy,y, w0,dw,w, vel,v, minx;
    float ang1, ang2, ang3, ang4, v1, v2, v3, v4, factor, *ctrace;
    bool angle, pass, cmplx;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    cmplx = (bool) (SF_COMPLEX == sf_gettype(in));

    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o1",&w0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("dim",&dim)) dim=2;
    /* [2,3] Dimensionality: filter 2-D planes or 3-D cubes */
    if (3==dim) {
	if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");
	if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
	if (!sf_histfloat(in,"o3",&y0)) sf_error("No o3= in input");
	n3 = sf_leftsize(in,3);		
    } else {
	ny=1;
	y0=0.;
	dy=dx;
	n3 = sf_leftsize(in,2);
    }
    ctrace = sf_floatalloc(cmplx? 2*nw: nw);

    if (!sf_getbool("angle",&angle)) angle=false;
    /* Filter based on angle (or velocity) */

    if (angle) {	
	if (!sf_getfloat("v",&v)) v=-1.;
	/* constant velocity (if angle-y)
	   The default is d(frequency)/d(wavenumber) */
	if (v < 0.) v=dw/hypotf(dx,dy);

	if (!sf_getfloat("ang1",&ang1)) ang1=-50.;
	if (!sf_getfloat("ang2",&ang2)) ang2=-45.;
	if (!sf_getfloat("ang3",&ang3)) ang3=45.;
	if (!sf_getfloat("ang4",&ang4)) ang4=50.;
	/* Angle gate (in degrees) */
	v1=tanf(SF_PI*ang1/180.)*v;
	v2=tanf(SF_PI*ang2/180.)*v;
	v3=tanf(SF_PI*ang3/180.)*v;
	v4=tanf(SF_PI*ang4/180.)*v;
    } else {
	if (!sf_getfloat("v1",&v1)) v1=0.;
	if (!sf_getfloat("v2",&v2)) v2=0.1;
	if (!sf_getfloat("v3",&v3)) v3=99999.;
	if (!sf_getfloat("v4",&v4)) v4=999999.;
	/* Velocity gate */
    }

    if ((v1>=v2)||(v2>=v3)||(v3>=v4)) 
	sf_error("Need v1 < v2 < v3 < v4, got %g ? %g ? %g ? %g",v1,v2,v3,v4); 

    if (!sf_getbool("pass",&pass)) pass=true;
    /* Pass or reject band */

    minx = hypotf(dx,dy)*FLT_EPSILON;

    for (i3=0; i3 < n3; i3++) { 
	for (iy=0; iy < ny; iy++) {
	    y = y0 + iy*dy;
	    for (ix=0; ix < nx; ix++) {
		x = x0 + ix*dx;
		x = hypotf(x,y)+minx;

		sf_floatread(ctrace,cmplx? 2*nw: nw,in);

		for (iw=0; iw < nw; iw++) {
		    w = w0+iw*dw;	    
		    vel = w/x; 

		    if ((vel>=v1) && (vel<=v2)) {
			factor=1.-sinf(0.5*SF_PI*(vel-v1)/(v2-v1));
		    } else if ((vel>=v2) && (vel<=v3)) {
			factor=0.;
		    } else if ((vel>=v3) && (vel<=v4)) {
			factor=sinf(0.5*SF_PI*(vel-v3)/(v4-v3));
		    } else {
			factor=1.;
		    }
			
		    if (pass) factor=1.-factor;
		    
		    if (cmplx) {
			ctrace[2*iw] *= factor;
			ctrace[2*iw+1] *= factor;
		    } else {
			ctrace[iw] *= factor;
		    }
		} /* iw */

		sf_floatwrite(ctrace,cmplx? 2*nw: nw,out);
	    } /* ix */
	} /* iy */
    } /* i3 */


    exit(0);
}

/* 	$Id$	 */
