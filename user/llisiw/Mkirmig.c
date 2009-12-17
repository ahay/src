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

static void doubint(int nt, float *trace)
/* causal and anticausal integratio */
{
    int it;
    float tt;

    tt = trace[0];
    for (it=1; it < nt; it++) {
	tt += trace[it];
	trace[it] = tt;
    }
    tt = trace[nt-1];
    for (it=nt-2; it >=0; it--) {
	tt += trace[it];
	trace[it] = tt;
    }
}
    
static float pick(float ti, float deltat, 
		  const float *trace,
		  int nt, float dt, float t0)
/* pick a traveltime sample from a trace */
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return 0.0;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return 0.0;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return 0.0;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    return imp*(
	2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	(1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	(1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
}

int main(int argc, char* argv[])
{
    int nt, nx, ny, ns, nh, nz, nzx, ix, ih, i, is, ist, iht;
    float *trace, *out, **table, *stable, *rtable;
    float ds, s0, x0, y0, dy, s, h,h0,dh,dx,ti,t0,t1,t2,dt,z0,dz,aal;
    sf_file inp, mig, tbl;

    sf_init (argc,argv);
    inp = sf_input("in");
    tbl = sf_input("table"); /* traveltime table */
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nh)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&ns)) sf_error("No n3=");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1=");

    if (!sf_histfloat(inp,"o2",&h0)) sf_error("No o2=");
    if (!sf_histfloat(inp,"d2",&dh)) sf_error("No d2=");

    if (!sf_histfloat(inp,"o3",&s0)) sf_error("No o3=");
    if (!sf_histfloat(inp,"d3",&ds)) sf_error("No d3=");

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
    sf_putfloat(mig,"o2",x0);
    sf_putfloat(mig,"d2",dx);

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */

    nzx = nz*nx;

    table = sf_floatalloc2(nzx,ny);
    sf_floatread(table[0],nzx*ny,tbl);

    out = sf_floatalloc(nzx);
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

	for (i=0; i < nzx; i++) {
	    out[i] = 0.;
	}
	
	for (ih=0; ih < nh; ih++) { /* offset */
	    h = h0+ih*dh;

	    /* place in the table */
	    /* nearest neighbor interpolation */

	    iht = 0.5 + (s+h-y0)/dy;
	    if (iht < 0) iht=0;
	    if (iht >= ny) iht=ny-1;

	    rtable = table[iht];

	    sf_floatread (trace,nt,inp);
	    doubint(nt,trace);

	    /* Add aperture limitation later */

	    for (ix=0; ix < nzx; ix++) { /* image */
		t1 = stable[ix];
		t2 = rtable[ix];
		ti = t1+t2;

		/* Add antialiasing later */
		out[ix] += pick(ti,0.,trace,nt,dt,t0);
	    } 
	} 
	
	sf_floatwrite(out,nzx,mig);        
    }

    exit(0);
}
