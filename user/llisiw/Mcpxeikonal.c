/* Iterative Complex Eikonal Solver */
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

#include "cpxeiko.h"
#include "fastmarchcpx.h"

int main(int argc, char* argv[])
{
    bool *k, velocity, verb;
    int dim, i, n[SF_MAX_DIM], *m, it, nt, iter, niter, cgiter, istep, nstep;
    float d[SF_MAX_DIM], o[SF_MAX_DIM], *s, *tr, *ti, *tr0, *ti0;
    float *w, *wr, *wi, *dw, *rhs, *x0, *fix, *exact;
    sf_complex *t, *t0;
    float rhsnorm, rhsnorm0, rhsnorm1, step;
    char key[4];
    sf_file in, out, vel, mask, ref, witer, dwiter, rhsiter;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    /* read input dimension */
    dim = sf_filedims(in,n);
    
    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; d[2] = d[1]; o[2] = o[1];
    }

    /* read initial guess */
    t  = sf_complexalloc(nt);
    tr = sf_floatalloc(nt);
    ti = sf_floatalloc(nt);
    
    sf_complexread(t,nt,in);
    
    for (it=0; it < nt; it++) {
	tr[it] = crealf(t[it]);
	ti[it] = cimagf(t[it]);
    }

    /* read background velocity */
    if (NULL == sf_getstring("vel"))
	sf_error("Need background velocity vel=");
    vel = sf_input("vel");

    s = sf_floatalloc(nt);
    sf_floatread(s,nt,vel);
    sf_fileclose(vel);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./s[it]*1./s[it];
	}
    }

    if(!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
    /* number of conjugate gradient iterations */
    
    if (!sf_getint("nstep",&nstep)) nstep=10;
    /* number of linesearch */
    
    /* output w at each iteration */
    if (NULL != sf_getstring("witer")) {
	witer = sf_output("witer");
	sf_settype(witer,SF_FLOAT);
	sf_putint(witer,"n3",n[2]);
	sf_putint(witer,"n4",niter+1);
    } else {
	witer = NULL;
    }

    /* output dw at each iteration */
    if (NULL != sf_getstring("dwiter")) {
	dwiter = sf_output("dwiter");
	sf_settype(dwiter,SF_FLOAT);
	sf_putint(dwiter,"n3",n[2]);
	sf_putint(dwiter,"n4",niter);
    } else {
	dwiter = NULL;
    }

    /* output rhs at each iteration */
    if (NULL != sf_getstring("rhsiter")) {
	rhsiter = sf_output("rhsiter");
	sf_settype(rhsiter,SF_FLOAT);
	sf_putint(rhsiter,"n3",n[2]);
	sf_putint(rhsiter,"n4",niter+1);
    } else {
	rhsiter = NULL;
    }
    
    /* read boundary condition */
    k = sf_boolalloc(nt);
    m = sf_intalloc(nt);

    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	sf_intread(m,nt,mask);
	sf_fileclose(mask);

	for (it=0; it < nt; it++) {
	    if (m[it] != 0)
		k[it] = true;
	    else
		k[it] = false;
	}
    } else {
	for (it=0; it < nt; it++) {
	    k[it] = false;
	}
    }

    t0    = sf_complexalloc(nt);
    tr0   = sf_floatalloc(n[1]);
    ti0   = sf_floatalloc(n[1]);
    fix   = sf_floatalloc(nt);
    exact = sf_floatalloc(nt);
    
    ref = sf_input("ref");
    sf_complexread(t0,nt,ref);
    sf_fileclose(ref);
    
    for (it=0; it < n[1]; it++) {
	tr0[it] = crealf(t0[it*n[0]]);
	ti0[it] = cimagf(t0[it*n[0]]);
    }
    
    for (it=0; it < nt; it++) {
	exact[it] = cimagf(t0[it]);
    }
    cpxeiko_ref(dim,n,d,exact,fix);

    free(m);
    free(t0);
    free(exact);
    m = NULL;
    t0 = NULL;
    exact = NULL;

    /* allocate temporary memory */
    w     = sf_floatalloc(nt);
    wr    = sf_floatalloc(nt);
    wi    = sf_floatalloc(nt);
    x0    = sf_floatalloc(nt);
    dw    = sf_floatalloc(nt);
    rhs   = sf_floatalloc(nt);

    /* initialize fastmarchcpx */
    fastmarchcpx_init(n,d);

    /* initialize cpxeiko */
    cpxeiko_init(dim,n,nt,d);

    /* initial misfit */
    cpxeiko_set(tr,ti);

    cpxeiko_forw(false,ti,w);
    for (it=0; it < nt; it++) x0[it] = fix[it]-w[it];
    
    cpxeiko_forw(false,tr,wr);
    cpxeiko_forw(true, ti,wi);
    for (it=0; it < nt; it++) rhs[it] = -wr[it]-wi[it];
    
    if (NULL != witer) sf_floatwrite(w,nt,witer);
    if (NULL != rhsiter) sf_floatwrite(rhs,nt,rhsiter);
    
    /* calculate L2 data-misfit */
    rhsnorm0 = cblas_snrm2(nt,rhs,1);
    rhsnorm1 = rhsnorm0;
    rhsnorm = rhsnorm0;
    
    sf_warning("Iteration 0 of %d:\t misfit=%g, w=%g.",niter,rhsnorm/rhsnorm0,cblas_snrm2(nt,w,1));

    for (iter=0; iter < niter; iter++) {

	/* clean-up */
	for (it=0; it < nt; it++) {
	    dw[it] = 0.;
	}
	
	/* solve dw */
	sf_solver(cpxeiko_operator,sf_cgstep,nt,nt,dw,rhs,cgiter,"known",k,"x0",x0,"verb",verb,"end");
	sf_cgstep_close();
	
	if (NULL != dwiter) sf_floatwrite(dw,nt,dwiter);

	/* linesearch */
	for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {

	    /* update real and imaginary slowness */
	    for (it=0; it < nt; it++) {
		wi[it] = w[it]+step*dw[it];
		if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
		
		wr[it] = s[it]+wi[it];
	    }
	    
	    /* forward fast-marching for stencil time */
	    fastmarchcpx(tr,tr0,wr);
	    fastmarchcpx(ti,ti0,wi);
	    
	    cpxeiko_set(tr,ti);
	    
	    cpxeiko_forw(false,tr,wr);
	    cpxeiko_forw(true, ti,wi);	    
	    for (it=0; it < nt; it++) rhs[it] = -wr[it]-wi[it];
	    
	    rhsnorm = cblas_snrm2(nt,rhs,1);
	    
	    /* break */
	    if (rhsnorm < rhsnorm1) {
		for (it=0; it < nt; it++) w[it] = wi[it];		
		rhsnorm1 = rhsnorm;
		break;
	    }
	}	
	
	cpxeiko_forw(false,ti,w);
	for (it=0; it < nt; it++) x0[it] = fix[it]-w[it];
	
	if (NULL != witer) sf_floatwrite(w,nt,witer);
	if (NULL != rhsiter) sf_floatwrite(rhs,nt,rhsiter);
	
	sf_warning("Iteration %d of %d:\t misfit=%g, w=%g, dw=%g (line-search %d).",iter+1,niter,rhsnorm/rhsnorm0,cblas_snrm2(nt,w,1),cblas_snrm2(nt,dw,1),istep);
    }

    for (it=0; it < nt; it++)
	t[it] = sf_cmplx(tr[it],ti[it]);

    sf_complexwrite(t,nt,out);
  
    exit(0);
}
