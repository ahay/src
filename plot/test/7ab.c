/* */
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
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    const int nt=300, nx=64, nz=300, nb=50; 
    int it, ib, iz, ix;
    float v=1.,t,z,x,x0, theta, *b, **tdat, **zdat;
    float top=3.2, c1=.9, c2=6.8, d1=0.02, d2=0.12, r;
    sf_file c, d;

    sf_init(argc,argv);
    
    b = sf_floatalloc(nb);
    tdat = sf_floatalloc2(nt,nx);
    zdat = sf_floatalloc2(nz,nx);

    if (!sf_getfloat("top",&top)) top=5.;
    if (!sf_getfloat("c1",&c1)) c1=0.5;
    if (!sf_getfloat("c2",&c2)) c2=5.;

    vp_init();

    vp_uorig (-c1,-.5);
    vp_uclip (0.,top-3.,4.,top);
    vp_umove (0.,top-0.);  
    vp_udraw (0.,top-4.);
    vp_umove (0.,top-0.);  
    vp_udraw (4.,top-0.);

    for (z=.4; z < 4.; z += .4) {
	vp_penup ();
	x0 = z * tanf( SF_PI*45./180.);
	for (x=0.; x < 4.; x += .01) {
	    t = hypotf(z,x-x0)/v;
	    vp_upendn (x,top-t);
	}
    }

    	
    for (ib=0; ib < nb; ib++) {
	b[ib] = expf(-3.*(ib+1.)/20.) * sinf(SF_PI*(ib+1.)/10.);
    }
    
    if (NULL != sf_getstring("c")) {
	c = sf_output("c");
	sf_setformat(c,"native_float");
	sf_putint(c,"n1",nt);
	sf_putint(c,"n2",nx);
	sf_putfloat(c,"d1",d1);
	sf_putfloat(c,"d2",d2);
	
	for (ix=0; ix < nx; ix++) {
	    for (it=0; it < nt; it++) {
		tdat[ix][it] = 0.;
	    }
	}
	
	for (iz=0; iz < 12; iz++) {
	    z = (iz+1.)*nt/12.;
	    x0 = z * tanf( SF_PI*45./180.);
	    for (ix=0; ix < nx; ix++) {
		x = (ix+1.) * d2/d1;
		t = hypotf(z,x-x0)/v;
		for (ib=0; ib < nb; ib++) {
		    it = t+ib-1;
		    if (it < nt) tdat[ix][it] += b[ib];
		}
	    }
	}

	sf_floatwrite(tdat[0],nt*nx,c);
    }

    vp_uorig (-c2,-.5);
    vp_uclip (0.,top-3.,4.,top);
    vp_umove (0.,top-0.);  
    vp_udraw (0.,top-4.);
    vp_umove (0.,top-0.);   
    vp_udraw (4.,top-0.);

    for(t=.4; t<6.; t += .4) {
	vp_penup ();
	x0 = t / sinf( SF_PI*45./180.);
	for (theta=-89.5; theta<89.5; theta += 1.) {
	    z =      t * cosf (SF_PI*theta/180.);
	    x = x0 + t * sinf (SF_PI*theta/180.);
	    vp_upendn (x,top-z);
	}
    }

    if (NULL != sf_getstring("d")) {
	d = sf_output("d");
	sf_setformat(d,"native_float");
	sf_putint(d,"n1",nz);
	sf_putint(d,"n2",nx);
	sf_putfloat(d,"d1",d1);
	sf_putfloat(d,"d1",d2);
	
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
	    zdat[ix][iz] = 0.;
	    }
	}
	
	for (it=0; it < 20; it++) {
	    t = (it+1.)*nz/20.;
	    x0 = t / sinf( SF_PI*45./180.);
	    for (theta=-89.5; theta<89.5; theta += 1.) {
		z =      t * cosf (SF_PI*theta/180.);
		x = x0 + t * sinf (SF_PI*theta/180.);
		ix = x * d1/d2;
		r = hypotf(z,x-x0);
		for (ib=0; ib < nb; ib++) {
		    iz = z + ib-1;
		    if (iz >= 0 && iz < nz && ix >=0 && ix < nx)
			zdat[ix][iz] += b[ib]*r;
		}
	    }
	}
	
	sf_floatwrite(zdat[0],nz*nx,d);
    }

    exit(0);
}

/* 	$Id: 7ab.c 7107 2011-04-10 02:04:14Z ivlad $	 */
