/* 2-D Prestack Kirchhoff depth migration. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "kirmig.h"
#include "tinterp.h"

int main(int argc, char* argv[])
{
    char *unit, *what, *type;
    bool cig;
    int nt, nx, ny, ns, nh, nz, nzx, ix, iz, ih, is, ist, iht;
    float *trace, **out, **table, *stable, *rtable, **tablex, *stablex, *rtablex;
    float ds, s0, x0, y0, dy, s, h, h0, dh, dx, ti, t0, t1, t2, dt, z0, dz, aal, tx, aper, thres;
    sf_file inp, mig, tbl, der;

    sf_init (argc,argv);
    inp = sf_input("in");
    tbl = sf_input("table"); /* traveltime table */
    der = sf_input("deriv"); /* source derivative table */
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&ns)) sf_error("No n3= in input");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(inp,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o3",&s0)) sf_error("No o3= in input");
    if (!sf_histfloat(inp,"d3",&ds)) sf_error("No d3= in input");

    if (1==nh) dh=0.0;
    if (1==ns) ds=0.0;

    if (!sf_histint(tbl,"n1",&nz)) sf_error("No n1= in table");
    if (!sf_histint(tbl,"n2",&nx)) sf_error("No n2= in table");
    if (!sf_histint(tbl,"n3",&ny)) sf_error("No n3= in table");

    if (!sf_histfloat(tbl,"o1",&z0)) sf_error("No o1= in table");
    if (!sf_histfloat(tbl,"d1",&dz)) sf_error("No d1= in table");

    if (!sf_histfloat(tbl,"o2",&x0)) sf_error("No o2= in table");
    if (!sf_histfloat(tbl,"d2",&dx)) sf_error("No d2= in table");

    if (!sf_histfloat(tbl,"o3",&y0)) sf_error("No o3= in table");
    if (!sf_histfloat(tbl,"d3",&dy)) sf_error("No d3= in table");

    if (!sf_getfloat("aperture",&aper)) aper=90.;
    /* migration aperture (in degree) */

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */
    
    if (!sf_getfloat("threshold",&thres)) thres=0.;
    /* source region thresholding */

    if (!sf_getbool("cig",&cig)) cig=false;
    /* y - output common offset gathers */

    sf_putint(mig,"n1",nz);
    sf_putint(mig,"n2",nx);

    sf_putfloat(mig,"o1",z0);
    sf_putfloat(mig,"d1",dz);
    sf_putfloat(mig,"o2",x0);
    sf_putfloat(mig,"d2",dx);
    sf_putstring(mig,"label1","Depth");
    sf_putstring(mig,"label2","Lateral");
    unit = sf_histstring(inp,"unit2");
    if (NULL != unit) sf_putstring(mig,"unit1",unit);
    if (cig) {
        sf_putint(mig,"n3",nh);
        sf_putfloat(mig,"o3",h0);
        sf_putfloat(mig,"d3",dh);
        sf_putstring(mig,"label3","Offset");
        if (NULL != unit) sf_putstring(mig,"unit3",unit);
    }

    nzx = nz*nx;

    /* read traveltime table */
    table = sf_floatalloc2(nzx,ny);
    sf_floatread(table[0],(off_t)nzx*(off_t)ny,tbl);
    sf_fileclose(tbl);

    /* read derivative table */
    tablex = sf_floatalloc2(nzx,ny);
    sf_floatread(tablex[0],(off_t)nzx*(off_t)ny,der);
    sf_fileclose(der);    

    out = sf_floatalloc2(nzx,cig ? nh : 1);
    trace = sf_floatalloc(nt);

    stable  = sf_floatalloc(nzx);
    stablex = sf_floatalloc(nzx);
    rtable  = sf_floatalloc(nzx);
    rtablex = sf_floatalloc(nzx);

    if (NULL == (type = sf_getstring("type"))) type="hermit";
    /* type of interpolation (default Hermit) */    

    if (NULL == (what = sf_getstring("what"))) what="expanded";
    /* Hermite basis functions (default expanded) */

    /* initialize interpolation */
    tinterp_init(nzx,dy,what);

    /* initialize summation */
    kirmig_init(nt,dt,t0);

    for (is=0; is < ns; is++) { /* shot */
	s = s0 + is*ds;

	/* cubic Hermite spline interpolation */
	ist = (s-y0)/dy;
	if (ist <= 0) {
	    for (ix=0; ix < nzx; ix++) {
		stable[ix]  = table[0][ix];
		stablex[ix] = tablex[0][ix];
	    }
	} else if (ist >= ny-1) {
	    for (ix=0; ix < nzx; ix++) {
		stable[ix]  = table[ny-1][ix];
		stablex[ix] = tablex[ny-1][ix];
	    }
	} else {
	    switch (type[0]) {
		case 'l': /* linear */
		    tinterp_linear(stable, s-ist*dy-y0,table[ist], table[ist+1]);
		    tinterp_linear(stablex,s-ist*dy-y0,tablex[ist],tablex[ist+1]);
		    break;

		case 'h': /* hermit */
		    tinterp_hermite(stable, s-ist*dy-y0,table[ist],table[ist+1],tablex[ist],tablex[ist+1]);
		    dinterp_hermite(stablex,s-ist*dy-y0,table[ist],table[ist+1],tablex[ist],tablex[ist+1]);
		    break;
	    }
	}

        if (!cig || 0 == is)
	    memset (&out[0][0], 0, cig ? nzx*nh*sizeof(float) : nzx*sizeof(float));

	for (ih=0; ih < nh; ih++) { /* offset */
	    h = h0+ih*dh;

	    /* cubic Hermite spline interpolation */
	    iht = (s+h-y0)/dy;
	    if (iht <= 0) {
		for (ix=0; ix < nzx; ix++) {
		    rtable[ix]  = table[0][ix];
		    rtablex[ix] = tablex[0][ix];
		}
	    } else if (iht >= ny-1) {
		for (ix=0; ix < nzx; ix++) {
		    rtable[ix]  = table[ny-1][ix];
		    rtablex[ix] = tablex[ny-1][ix];
		}
	    } else {
		switch (type[0]) {
		case 'l': /* linear */
		    tinterp_linear(rtable, s+h-iht*dy-y0,table[iht], table[iht+1]);
		    tinterp_linear(rtablex,s+h-iht*dy-y0,tablex[iht],tablex[iht+1]);
		    break;

		case 'h': /* hermit */
		    tinterp_hermite(rtable, s+h-iht*dy-y0,table[iht],table[iht+1],tablex[iht],tablex[iht+1]);
		    dinterp_hermite(rtablex,s+h-iht*dy-y0,table[iht],table[iht+1],tablex[iht],tablex[iht+1]);
		}
	    }

	    /* read trace */
	    sf_floatread (trace,nt,inp);
	    doubint(nt,trace);

	    for (ix=0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) { /* image */
		    /* aperture (cone angle) */
		    if (h >= 0.) {
			if (atanf((s-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-s-h)/(iz*dz))*180./SF_PI > aper) continue;
		    } else {
			if (atanf((s+h-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-s)/(iz*dz))*180./SF_PI > aper) continue;
		    }

		    t1 = stable[ix*nz+iz];
		    t2 = rtable[ix*nz+iz];
		    ti = t1+t2;
		    
		    if (ti < thres) continue;

		    tx = SF_MAX(fabsf(stablex[ix*nz+iz]*ds),fabsf(rtablex[ix*nz+iz]*dh));
		    pick(true,ti,tx*aal,out[cig ? ih : 0]+ix*nz+iz,trace);
		} 
	    }
	}
        if (!cig) 
            sf_floatwrite(out[0],(off_t)nzx,mig);
    }
    if (cig)
        sf_floatwrite(out[0],(off_t)nzx*(off_t)nh,mig);

    exit(0);
}
