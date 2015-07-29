/* 2-D Prestack Kirchhoff depth migration. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "kkirmig.h"

int main(int argc, char* argv[])
{
    char *unit;
    bool cig;
    int nt, nx, ny, ns, nh, nz, nzx, ix, ih, is, ist, iht;
    float *trace, **out, **table, *stable, *rtable, **tablex, *stablex, *rtablex;
    float ds, s0, x0, y0, dy, s, h,h0,dh,dx,ti,t0,t1,t2,dt,z0,dz,aal, tx;
    sf_file inp, mig, tbl;

    sf_init (argc,argv);
    inp = sf_input("in");
    tbl = sf_input("table"); /* traveltime table */
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

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */
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

    if (NULL != sf_getstring("tablex")) {
	tbl = sf_input("tablex");

	/* table of dT/ds */
	tablex = sf_floatalloc2(nzx,ny);
	sf_floatread(tablex[0],(off_t)nzx*(off_t)ny,tbl);
	sf_fileclose(tbl);
    } else {
	tablex = NULL;
    }

    out = sf_floatalloc2(nzx,cig ? nh : 1);
    trace = sf_floatalloc(nt);

    for (is=0; is < ns; is++) { /* shot */
	s = s0 + is*ds;

	/* Add accurate interpolation later */

	/* place in the table */
	/* nearest neighbor interpolation */

	ist = 0.5 + (s-y0)/dy;
	if (ist < 0) ist=0;
	if (ist >= ny) ist=ny-1;

	stable = table[ist];
	stablex = (NULL==tablex)? NULL:tablex[ist];

        if (!cig || 0 == is)
	    memset (&out[0][0], 0, cig ? nzx*nh*sizeof(float) : nzx*sizeof(float));

	for (ih=0; ih < nh; ih++) { /* offset */
	    h = h0+ih*dh;

	    /* place in the table */
	    /* nearest neighbor interpolation */

	    iht = 0.5 + (s+h-y0)/dy;
	    if (iht < 0) iht=0;
	    if (iht >= ny) iht=ny-1;

	    rtable = table[iht];
	    rtablex = (NULL==tablex)? NULL:tablex[iht];

	    sf_floatread (trace,nt,inp);
	    doubint(nt,trace);

	    /* Add aperture limitation later */

	    for (ix=0; ix < nzx; ix++) { /* image */
		t1 = stable[ix];
		t2 = rtable[ix];
		ti = t1+t2;
		
		tx = (NULL==tablex)? 0.: SF_MAX(fabsf(stablex[ix]*ds),
						fabsf(rtablex[ix]*dh));
                out[cig ? ih : 0][ix] += pick(ti,tx*aal,trace,nt,dt,t0);
	    } 
	}
        if (!cig) 
            sf_floatwrite(out[0],(off_t)nzx,mig);
    }
    if (cig)
        sf_floatwrite(out[0],(off_t)nzx*(off_t)nh,mig);

    exit(0);
}

