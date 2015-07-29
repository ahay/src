/* Iterative solution of the linearized eikonal equation. */
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

#include "upgrad.h"

int main(int argc, char* argv[])
{
    bool adj, inv, squared;
    int dim, i, n[SF_MAX_DIM], it, nt, niter, iter, *m;
    float d[SF_MAX_DIM], *t, *dt, *s, *ds, *t0, tol;
    double err;
    char key[4];
    const char *what;
    upgrad upg;
    sf_file dtime, time, slow, mask, time0;
    
    sf_init(argc,argv);
    time = sf_input("in");
    dtime = sf_output("out");
    
    dim = sf_filedims(time,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(time,key,d+i)) sf_error("No %s= in input",key);
	nt *= n[i];
    }

    if (NULL == (what = sf_getstring("what"))) what="time";
    /* what to compute */

    if (!sf_getbool("squared",&squared)) squared=true;
    /* if slowness is squared */

    upg = upgrad_init(dim,n,d);

    t = sf_floatalloc(nt);
    s = sf_floatalloc(nt);

    sf_floatread(t,nt,time);

    switch (what[0]) {
	case 's': /* slowness squared */
	    upgrad_set(upg,t);
	    upgrad_forw(upg,t,s); /* s = grad(t)*grad(t) */

	    sf_floatwrite(s,nt,dtime);
	    break;
	case 'l': /* linear operator */
	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */
	
	    if (!sf_getbool("inv",&inv)) inv=false;
	    /* inverse flag (for what=linear) */

	    /* get t0 */
	    if (NULL == sf_getstring("time"))
		sf_error("Need time=");

	    time0 = sf_input("time");
	    t0 = sf_floatalloc(nt);
	    sf_floatread(t0,nt,time0);
	    sf_fileclose(time0);
	    
	    upgrad_set(upg,t0);
	    
	    if (inv) {
		if (adj) {
		    upgrad_inverse(upg,s,t,NULL);
		} else {
		    upgrad_solve(upg,t,s,NULL);
		}
	    } else {
		if (adj) {
		    upgrad_adj(upg,s,t);
		} else {
		    upgrad_forw(upg,t,s);
		}
	    }

	    sf_floatwrite(s,nt,dtime);
	    break;
	case 'd': /* differentiation */
	    time0 = sf_input("time");
	    t0 = sf_floatalloc(nt);
	    sf_floatread(t0,nt,time0);
	    sf_fileclose(time0);

	    upgrad_set(upg,t);
	    upgrad_forw(upg,t0,s); /* s = grad(t0)*grad(t) */

	    if (!squared) {
		ds = sf_floatalloc(nt);
		upgrad_forw(upg,t,ds);
		for (it=0; it < nt; it++) {
		    s[it] /= sqrtf(ds[it]);
		}
		free(ds);
	    }
	    
	    sf_floatwrite(s,nt,dtime);
	    break;
	case 'i': /* integration */
	    /* boundary conditions */
	    if (NULL != sf_getstring("time")) {
		time0 = sf_input("time");
		t0 = sf_floatalloc(nt);
		sf_floatread(t0,nt,time0);
		sf_fileclose(time0);
	    } else {
		t0 = NULL;
	    }

	    /* get s */
	    slow = sf_input("slow");  /* slowness squared */
	    sf_floatread(s,nt,slow);
	    sf_fileclose(slow);

	    upgrad_set(upg,t);

	    if (!squared) {
		ds = sf_floatalloc(nt);
		upgrad_forw(upg,t,ds);
		for (it=0; it < nt; it++) {
		    s[it] *= sqrtf(ds[it]);
		}
		free(ds);
	    }

	    upgrad_solve(upg,s,t,t0); /* solve for t */
	    
	    sf_floatwrite(t,nt,dtime);
	    break;
	case 't': /* time */
	    slow = sf_input("slow"); /* slowness */
	    
	    if (!sf_getint("niter",&niter)) niter=1;
	    /* maximum number of iterations */
	    
	    if (!sf_getfloat("tol",&tol)) tol=0.001;
	    /* tolerance for convergence */
	    
	    ds = sf_floatalloc(nt);
	    dt = sf_floatalloc(nt);
	    m = sf_intalloc(nt);
	    
	    if (NULL != sf_getstring("mask")) {
		mask = sf_input("mask");
		if (SF_INT != sf_gettype(mask)) sf_error("Need int mask");
		sf_intread(m,nt,mask);
		sf_fileclose(mask);
	    } else {
		for (it=0; it < nt; it++) {
		    m[it]=1;
		}
	    }
	    
	    sf_floatread(s,nt,slow);
	    sf_fileclose(slow);
	    
	    for (iter=0; iter < niter; iter++) {
		upgrad_set(upg,t);
		upgrad_forw(upg,t,ds);
		
		for (it=0; it < nt; it++) {
		    ds[it] = m[it]? sqrtf(ds[it])*s[it]-ds[it]:0.0;
		}
		
		upgrad_solve(upg,ds,dt,NULL);
		
		for (it=0; it < nt; it++) {
		    t[it] += dt[it];
		}
		
		err = sqrt(cblas_dsdot(nt,dt,1,dt,1)/nt);
		if (err < tol) break;
		
		sf_warning("cycle=%d dt=%g",iter+1,err);
	    }
	    
	    sf_floatwrite(t,nt,dtime);
	    break;
    }

    exit(0);
}

