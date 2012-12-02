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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kirmig.h"
#include "tinterp.h"

int main(int argc, char* argv[])
{
    char *unit, *type;
    bool adj, cig, cmp;
    off_t nzx;
    int nt, nx, ny, ns, nh, nz, i, ix, iz, ih, is, ist, iht, ng;
    float *trace, **out, **table, *stable, *rtable, **tablex, *stablex, *rtablex;
    float ds, s0, x0, y0, dy, s, h, h0, dh, dx, ti, t0, t1, t2, dt, z0, dz, tau;
    float aal, tx, aper;
    sf_file dat, mig, tbl, der;

    sf_init (argc,argv);

    if (!sf_getbool("adj",&adj)) adj=true;
    /* y for migration, n for modeling */

    if (!sf_getbool("cmp",&cmp)) cmp=true;
    /* y for CMP gather, n for shot gather */    

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
	if (!sf_histint(dat,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(dat,"n2",&nh)) sf_error("No n2= in input");
	if (!sf_histint(dat,"n3",&ns)) sf_error("No n3= in input");
	
	if (!sf_histfloat(dat,"o1",&t0)) sf_error("No o1= in input");
	if (!sf_histfloat(dat,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(dat,"o2",&h0)) sf_error("No o2= in input");
	if (!sf_histfloat(dat,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(dat,"o3",&s0)) sf_error("No o3= in input");
	if (!sf_histfloat(dat,"d3",&ds)) sf_error("No d3= in input");
    } else {	
	if (!sf_getint("nt",&nt)) sf_error("Need nt="); /* time samples */
	if (!sf_getint("nh",&nh)) nh=1; /* offset/receiver samples */
	if (!sf_getint("ns",&ns)) ns=1; /* shot samples */
	
	if (!sf_getfloat("t0",&t0)) t0=0.0; /* time origin */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt="); /* time sampling */
	if (!sf_getfloat("h0",&h0)) h0=0.0; /* offset/receiver origin */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh="); /* offset/receiver sampling */
	if (!sf_getfloat("s0",&s0)) s0=0.0; /* shot origin */
	if (!sf_getfloat("ds",&ds)) sf_error("Need ds="); /* shot sampling */
    }

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

    if (!sf_getfloat("tau",&tau)) tau=0.;
    /* static time-shift (in second) */

    if (!sf_getfloat("aperture",&aper)) aper=90.;
    /* migration aperture (in degree) */

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */
    
    if (!sf_getbool("cig",&cig)) cig=false;
    /* y - output common offset/receiver gathers */

    ng = cig? nh: 1;

    if (adj) {
	sf_putint(mig,"n1",nz);
	sf_putint(mig,"n2",nx);
	
	sf_putfloat(mig,"o1",z0);
	sf_putfloat(mig,"d1",dz);
	sf_putfloat(mig,"o2",x0);
	sf_putfloat(mig,"d2",dx);
	sf_putstring(mig,"label1","Depth");
	sf_putstring(mig,"label2","Lateral");
	unit = sf_histstring(dat,"unit2");
	if (NULL != unit) sf_putstring(mig,"unit1",unit);
	if (cig) {
	    sf_putint(mig,"n3",nh);
	    sf_putfloat(mig,"o3",h0);
	    sf_putfloat(mig,"d3",dh);

	    if (cmp)
		sf_putstring(mig,"label3","Offset");
	    else
		sf_putstring(mig,"label3","Receiver");

	    if (NULL != unit) sf_putstring(mig,"unit3",unit);
	} else {
	    sf_putint(mig,"n3",1);
	}
    } else {
	sf_putint(dat,"n1",nt);
	sf_putint(dat,"n2",nh);
	sf_putint(dat,"n3",ns);
	
	sf_putfloat(dat,"o1",t0);
	sf_putfloat(dat,"d1",dt);

	sf_putfloat(dat,"o2",h0);
	sf_putfloat(dat,"d2",dh);

	sf_putfloat(dat,"o3",s0);
	sf_putfloat(dat,"d3",ds);

	sf_putstring(dat,"label1","Time");
	sf_putstring(dat,"unit1","s");	

	if (cmp)
	    sf_putstring(dat,"label2","Offset");
	else
	    sf_putstring(dat,"label2","Receiver");
	sf_putstring(dat,"label3","Shot");
    }
	
    nzx = (off_t) nz* (off_t) nx;

    /* read traveltime table */
    table = sf_floatalloc2(nzx,ny);
    sf_floatread(table[0],nzx*ny,tbl);
    sf_fileclose(tbl);

    /* read derivative table */
    tablex = sf_floatalloc2(nzx,ny);
    sf_floatread(tablex[0],nzx*ny,der);
    sf_fileclose(der);    

    out = sf_floatalloc2(nzx,ng);
    trace = sf_floatalloc(nt);

    stable  = sf_floatalloc(nzx);
    stablex = sf_floatalloc(nzx);
    rtable  = sf_floatalloc(nzx);
    rtablex = sf_floatalloc(nzx);

    if (NULL == (type = sf_getstring("type"))) type="hermit";
    /* type of interpolation (default Hermit) */    

    /* initialize interpolation */
    tinterp_init(nzx,dy);

    /* initialize summation */
    kirmig_init(nt,dt,t0);
    
    if (adj) {
	memset (out[0], 0, nzx*ng*sizeof(float));
    } else {		
	sf_floatread(out[0],nzx,mig);
    }

    for (is=0; is < ns; is++) { /* shot */
	s = s0+is*ds;

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

	for (ih=0; ih < nh; ih++) { /* offset/receiver */
	    h = h0+ih*dh;

	    /* cubic Hermite spline interpolation */
	    iht = cmp? (s+h-y0)/dy: (h-y0)/dy;
	    if (iht <= 0) {
		for (i=0; i < nzx; i++) {
		    rtable[i]  = table[0][i];
		    rtablex[i] = tablex[0][i];
		}
	    } else if (iht >= ny-1) {
		for (i=0; i < nzx; i++) {
		    rtable[i]  = table[ny-1][i];
		    rtablex[i] = tablex[ny-1][i];
		}
	    } else {
		switch (type[0]) {
		    case 'l': /* linear */
			tinterp_linear(rtable, cmp? s+h-iht*dy-y0: h-iht*dy-y0,table[iht],table[iht+1]);
			dinterp_linear(rtablex,cmp? s+h-iht*dy-y0: h-iht*dy-y0,table[iht],table[iht+1]);
			break;
			
		    case 'p': /* partial */
			tinterp_partial(rtable, cmp? s+h-iht*dy-y0: h-iht*dy-y0,nz,nx,dx,table[iht],table[iht+1]);
			dinterp_partial(rtablex,cmp? s+h-iht*dy-y0: h-iht*dy-y0,nz,nx,dx,table[iht],table[iht+1]);
			break;

		    case 'h': /* hermit */
			tinterp_hermite(rtable, cmp? s+h-iht*dy-y0: h-iht*dy-y0,table[iht],table[iht+1],tablex[iht],tablex[iht+1]);
			dinterp_hermite(rtablex,cmp? s+h-iht*dy-y0: h-iht*dy-y0,table[iht],table[iht+1],tablex[iht],tablex[iht+1]);
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
#pragma omp parallel for private(iz,ix,t1,t2,ti,tx)
#endif
	    for (i=0; i < nzx; i++) { 
		iz = i%nz;
		ix = (i-iz)/nz;

		/* aperture (cone angle) */
		if (cmp) {
		    if (h >= 0.) {
			if (atanf((s-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-s-h)/(iz*dz))*180./SF_PI > aper) continue;
		    } else {
			if (atanf((s+h-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-s)/(iz*dz))*180./SF_PI > aper) continue;
		    }
		} else {
		    if (h-s >= 0.) {
			if (atanf((s-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-h)/(iz*dz))*180./SF_PI > aper) continue;
		    } else {
			if (atanf((h-x0-ix*dx)/(iz*dz))*180./SF_PI > aper) continue;
			if (atanf((x0+ix*dx-s)/(iz*dz))*180./SF_PI > aper) continue;
		    }
		}

		t1 = stable[i];
		t2 = rtable[i];
		ti = t1+t2+tau;

		tx = SF_MAX(fabsf(stablex[i]*ds),fabsf(rtablex[i]*dh));
		pick(adj,ti,tx*aal,out[cig ? ih : 0]+i,trace);
	    }

	    if (!adj) {
		doubint(nt,trace);
		sf_floatwrite (trace,nt,dat);
	    }
	} /* ih */
    }
    
    if (adj) sf_floatwrite(out[0],nzx*ng,mig);

    exit(0);
}
