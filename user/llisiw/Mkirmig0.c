/* 2-D Post-stack Kirchhoff depth migration. */
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kirmig.h"
#include "tinterp.h"

int main(int argc, char* argv[])
{
    bool adj;
    int nt, nx, ny, ns, nz, nzx, ix, iz, i, is, ist;
    float *trace, *out, **table, **tablex, *stable, *stablex;
    float ds, s0, x0, y0, dy, s, dx, ti, t0, dt, z0, dz, aal, tx, aper;
    char *unit, *type;
    sf_file dat, mig, tbl, der;

    sf_init (argc,argv);

    if (!sf_getbool("adj",&adj)) adj=true;
    /* y for migration, n for modeling */

    if (adj) {
	dat = sf_input("in");
	mig = sf_output("out");
    } else {
	mig = sf_input("in");
	dat = sf_output("out");
    }

    tbl = sf_input("table"); /* traveltime table */
    der = sf_input("deriv"); /* source derivative table */

    if (adj) {
	if (!sf_histint(dat,"n1",&nt)) sf_error("No n1=");
	if (!sf_histint(dat,"n2",&ns)) sf_error("No n2=");
	
	if (!sf_histfloat(dat,"o1",&t0)) sf_error("No o1=");
	if (!sf_histfloat(dat,"d1",&dt)) sf_error("No d1=");
	
	if (!sf_histfloat(dat,"o2",&s0)) sf_error("No o2=");
	if (!sf_histfloat(dat,"d2",&ds)) sf_error("No d2=");
    } else {
	if (!sf_getint("nt",&nt)) sf_error("Need nt="); /* time samples */
	if (!sf_getint("ns",&ns)) sf_error("Need ns="); /* midpoint samples */
	
	if (!sf_getfloat("t0",&t0)) t0=0.0;                /* time origin */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");  /* time sampling */
	
	if (!sf_getfloat("s0",&s0)) s0=0.0;                /* midpoint origin */
	if (!sf_getfloat("ds",&ds)) sf_error("Need ds=");  /* midpoint sampling */
    }

    if (!sf_histint(tbl,"n1",&nz)) sf_error("No n1= in table");
    if (!sf_histint(tbl,"n2",&nx)) sf_error("No n2= in table");
    if (!sf_histint(tbl,"n3",&ny)) sf_error("No n3= in table");

    if (!sf_histfloat(tbl,"o1",&z0)) sf_error("No o1= in table");
    if (!sf_histfloat(tbl,"d1",&dz)) sf_error("No d1= in table");
    if (!sf_histfloat(tbl,"o2",&x0)) sf_error("No o2= in table");
    if (!sf_histfloat(tbl,"d2",&dx)) sf_error("No d2= in table");
    if (!sf_histfloat(tbl,"o3",&y0)) sf_error("No o3= in table");
    if (!sf_histfloat(tbl,"d3",&dy)) sf_error("No d3= in table");

    if (adj) {
	sf_putint(mig,"n1",nz);
	sf_putint(mig,"n2",nx);

	sf_putfloat(mig,"o1",z0);
	sf_putfloat(mig,"d1",dz);
	sf_putstring(mig,"label1","Depth");
	unit = sf_histstring(dat,"unit2");
	if (NULL != unit) sf_putstring(mig,"unit1",unit);

	sf_putfloat(mig,"o2",x0);
	sf_putfloat(mig,"d2",dx);
	sf_putstring(mig,"label2","Distance");
    } else {
	sf_putint(dat,"n1",nt);
	sf_putint(dat,"n2",ns);
	
	sf_putfloat(dat,"o1",t0);
	sf_putfloat(dat,"d1",dt);

	sf_putfloat(dat,"o2",s0);
	sf_putfloat(dat,"d2",ds);

	sf_putstring(dat,"label1","Time");
	sf_putstring(dat,"unit1","s");	
    }

    if (!sf_getfloat("aperture",&aper)) aper=90.;
    /* migration aperture (in degree) */

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */

    nzx = nz*nx;

    /* read traveltime table */
    table = sf_floatalloc2(nzx,ny);
    sf_floatread(table[0],nzx*ny,tbl);
    sf_fileclose(tbl);

    /* read derivative table */
    tablex = sf_floatalloc2(nzx,ny);
    sf_floatread(tablex[0],nzx*ny,der);
    sf_fileclose(der);

    out = sf_floatalloc(nzx);
    trace = sf_floatalloc(nt);

    stable  = sf_floatalloc(nzx);
    stablex = sf_floatalloc(nzx);

    if (NULL == (type = sf_getstring("type"))) type="hermit";
    /* type of interpolation (default Hermit) */

    /* initialize interpolation */
    tinterp_init(nzx,dy);

    if (adj) {
	for (i=0; i < nzx; i++) {
	    out[i] = 0.;
	}
    } else {
	sf_floatread(out,nzx,mig);
    }       

    kirmig_init(nt,dt,t0);

    for (is=0; is < ns; is++) { /* surface location */
	s = s0 + is*ds;

	/* cubic Hermite spline interpolation */
	ist = (s-y0)/dy;
	if (ist <= 0) {
	    for (i=0; i < nzx; i++) {
		stable[i]  = table[0][i];
		stablex[i] = tablex[0][i];
	    }
	} else if (ist >= ny-1) {
	    for (i=0; i < nzx; i++) {
		stable[i]  = table[ny-1][i];
		stablex[i] = tablex[ny-1][i];
	    }
	} else {
	    switch (type[0]) {
		case 'l': /* linear */
		    tinterp_linear(stable, s-ist*dy-y0,table[ist],table[ist+1]);
		    dinterp_linear(stablex,s-ist*dy-y0,table[ist],table[ist+1]);
		    break;

		case 'p': /* partial */
		    tinterp_partial(stable, s-ist*dy-y0,nz,nx,dx,table[ist],table[ist+1]);
		    dinterp_partial(stablex,s-ist*dy-y0,nz,nx,dx,table[ist],table[ist+1]);
		    break;

		case 'h': /* hermit */
		    tinterp_hermite(stable, s-ist*dy-y0,table[ist],table[ist+1],tablex[ist],tablex[ist+1]);
		    dinterp_hermite(stablex,s-ist*dy-y0,table[ist],table[ist+1],tablex[ist],tablex[ist+1]);
		    break;
	    }
	}	

	if (adj) {
	    /* read trace */
	    sf_floatread (trace,nt,dat);
	    doubint(nt,trace);
	} else {
	    for (i=0; i < nt; i++) {
		trace[i]=0.;
	    }
	}

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,ti,tx)
#endif
	for (i=0; i < nzx; i++) { 
	    iz = i%nz;
	    ix = (i-iz)/nz;

	    /* aperture (cone angle) */
	    if (fabsf(atanf((x0+ix*dx-s)/(iz*dz)))*180./SF_PI > aper)
		continue;

	    ti = 2.*stable[i];
	    tx = 2.*stablex[i];
	    
	    pick(adj,ti,fabsf(tx*ds*aal),out+i,trace);
	}
	
	if (!adj) {
	    doubint(nt,trace);
	    sf_floatwrite (trace,nt,dat);
	}
    }

    if (adj) sf_floatwrite(out,nzx,mig);        

    exit(0);
}
