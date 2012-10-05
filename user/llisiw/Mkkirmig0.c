/* 2-D Post-stack Kirchhoff depth migration. */
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
    int nt, nx, ny, ns, nz, nzx, ix, i, is, ist;
    float *trace, *out, **table, **tablex, *stable, *stablex;
    float ds, s0, x0, y0, dy, s, dx,ti,t0,dt,z0,dz,aal, tx;
    char *unit;
    sf_file inp, mig, tbl;

    sf_init (argc,argv);
    inp = sf_input("in");
    tbl = sf_input("table"); /* traveltime table */
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&ns)) sf_error("No n2=");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1=");

    if (!sf_histfloat(inp,"o2",&s0)) sf_error("No o2=");
    if (!sf_histfloat(inp,"d2",&ds)) sf_error("No d2=");

    if (!sf_histint(tbl,"n1",&nz)) sf_error("No n1= in table");
    if (!sf_histint(tbl,"n2",&nx)) sf_error("No n2= in table");
    if (!sf_histint(tbl,"n3",&ny)) sf_error("No n3= in table");

    if (!sf_histfloat(tbl,"o1",&z0)) sf_error("No o1= in table");
    if (!sf_histfloat(tbl,"d1",&dz)) sf_error("No d1= in table");
    if (!sf_histfloat(tbl,"o2",&x0)) sf_error("No o2= in table");
    if (!sf_histfloat(tbl,"d2",&dx)) sf_error("No d2= in table");
    if (!sf_histfloat(tbl,"o3",&y0)) sf_error("No o3= in table");
    if (!sf_histfloat(tbl,"d3",&dy)) sf_error("No d3= in table");

    sf_putint(mig,"n1",nz);
    sf_putint(mig,"n2",nx);

    sf_putfloat(mig,"o1",z0);
    sf_putfloat(mig,"d1",dz);
    sf_putstring(mig,"label1","Depth");
    unit = sf_histstring(inp,"unit2");
    if (NULL != unit) sf_putstring(mig,"unit1",unit);

    sf_putfloat(mig,"o2",x0);
    sf_putfloat(mig,"d2",dx);
    sf_putstring(mig,"label2","Distance");

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */

    nzx = nz*nx;

    table = sf_floatalloc2(nzx,ny);
    sf_floatread(table[0],nzx*ny,tbl);
    sf_fileclose(tbl);

    if (NULL != sf_getstring("tablex")) {
	tbl = sf_input("tablex");

	tablex = sf_floatalloc2(nzx,ny);
	sf_floatread(tablex[0],nzx*ny,tbl);
	sf_fileclose(tbl);
    } else {
	tablex = NULL;
    }

    out = sf_floatalloc(nzx);
    trace = sf_floatalloc(nt);

    for (i=0; i < nzx; i++) {
	out[i] = 0.;
    }

    for (is=0; is < ns; is++) { /* surface location */
	s = s0 + is*ds;

	/* Add accurate interpolation later */

	/* place in the table */
	/* nearest neighbor interpolation */

	ist = 0.5 + (s-y0)/dy;
	if (ist < 0) ist=0;
	if (ist >= ny) ist=ny-1;

	stable = table[ist];
	stablex = (NULL==tablex)? NULL:tablex[ist];

	sf_floatread (trace,nt,inp);
	doubint(nt,trace);

	/* Add aperture limitation later */

	for (ix=0; ix < nzx; ix++) { /* image */
	    ti = 2*stable[ix];
	    tx = (NULL==stablex)?0.:2*stablex[ix];
	    
	    out[ix] += pick(ti,fabsf(tx*ds*aal),trace,nt,dt,t0);
	}
    }

    sf_floatwrite(out,nzx,mig);        

    exit(0);
}
