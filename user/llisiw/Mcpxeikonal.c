/* Iterative complex eikonal solver */
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
    bool *known, velocity, verb;
    int dim, i, n[SF_MAX_DIM], *m, it, nt, iter, niter, cgiter, istep, nstep, count;
    float d[SF_MAX_DIM], o[SF_MAX_DIM], *s, *tr, *ti, *tr0, *ti0;
    float *w, *wr, *wi, *dw, *rhs, *rhsr, *rhsi, *x0, *fix, *wt, gama, ratio, tol, ***oper, *unit;
    sf_complex *t, *t0, **pdir;
    float rhsnorm, rhsnorm0, rhsnorm1, step;
    char key[4];
    sf_file in, out, vel, mask, ref, gap;
    sf_file witer, dwiter, rhsiter, upiter, x0iter, titer, wtiter, liniter, operiter;
    
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
    
    /* output stencil at each iteration */
    /* NOTE: only first iteration */
    if (NULL != sf_getstring("upiter")) {
	upiter = sf_output("upiter");
	sf_putint(upiter,"n3",n[2]);
	sf_putint(upiter,"n4",dim);
	sf_putint(upiter,"n5",1);
	pdir = sf_complexalloc2(nt,dim);
    } else {
	upiter = NULL;
	pdir = NULL;
    }
    
    /* output operator at each iteration */
    /* NOTE: only first iteration */
    if (NULL != sf_getstring("operiter")) {
	operiter = sf_output("operiter");
	sf_settype(operiter,SF_FLOAT);
	sf_putint(operiter,"n1",nt);
	sf_putint(operiter,"n2",nt);
	sf_putint(operiter,"n3",2);
	sf_putint(operiter,"n4",1);
	oper = sf_floatalloc3(nt,nt,2);
	unit = sf_floatalloc(nt);
    } else {
	operiter = NULL;
	oper = NULL;
	unit = NULL;
    }
    
    /* output x0 at each iteration */
    if (NULL != sf_getstring("x0iter")) {
	x0iter = sf_output("x0iter");
	sf_settype(x0iter,SF_FLOAT);
	sf_putint(x0iter,"n3",n[2]);
	sf_putint(x0iter,"n4",niter+1);
    } else {
	x0iter = NULL;
    }

    /* output linearized cost function at each iteration */
    if (NULL != sf_getstring("liniter")) {
	liniter = sf_output("liniter");
	sf_settype(liniter,SF_FLOAT);
	sf_putint(liniter,"n3",n[2]);
	sf_putint(liniter,"n4",niter);
    } else {
	liniter = NULL;
    }

    /* output t at each iteration */
    if (NULL != sf_getstring("titer")) {
	titer = sf_output("titer");
	sf_putint(titer,"n3",n[2]);
	sf_putint(titer,"n4",niter+1);
    } else {
	titer = NULL;
    }

    /* output wt at each iteration */
    if (NULL != sf_getstring("wtiter")) {
	wtiter = sf_output("wtiter");
	sf_putint(wtiter,"n3",n[2]);
	sf_putint(wtiter,"n4",niter);
    } else {
	wtiter = NULL;
    }

    /* read boundary condition */
    m = sf_intalloc(nt);
    known = sf_boolalloc(nt);

    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	sf_intread(m,nt,mask);
	sf_fileclose(mask);
	
	for (it=0; it < nt; it++) {
	    if (m[it] != 0)
		known[it] = true;
	    else
		known[it] = false;
	}
    } else {
	for (it=0; it < nt; it++) {
	    known[it] = false;
	}
    }

    t0  = sf_complexalloc(nt);
    tr0 = sf_floatalloc(nt);
    ti0 = sf_floatalloc(nt);
    fix = sf_floatalloc(nt);
    
    ref = sf_input("ref");
    sf_complexread(t0,nt,ref);
    sf_fileclose(ref);
    
    for (it=0; it < nt; it++) {
	tr0[it] = crealf(t0[it]);
	ti0[it] = cimagf(t0[it]);
    }
    
    cpxeiko_ref(dim,n,d,ti0,fix);

    free(m);
    free(t0);
    m = NULL;
    t0 = NULL;

    /* read model weighting */
    wt = sf_floatalloc(nt);

    if (NULL != sf_getstring("gap")) {
	gap = sf_input("gap");
	sf_floatread(wt,nt,gap);
	sf_fileclose(gap);
    } else {
	for (it=0; it < nt; it++) {
	    wt[it] = 1.;
	}
    }
    
    /* allocate temporary memory */
    w     = sf_floatalloc(nt);
    rhsr    = sf_floatalloc(nt);
    rhsi    = sf_floatalloc(nt);
    wr   = sf_floatalloc(nt);
    wi   = sf_floatalloc(nt);
    x0    = sf_floatalloc(nt);
    dw    = sf_floatalloc(nt);
    rhs   = sf_floatalloc(nt);

    /* initialize fastmarchcpx */
    /* NOTE: default accuracy 2nd order */
    fastmarchcpx_init(n,o,d,2);

    /* initialize cpxeiko */
    cpxeiko_init(dim,n,nt,d);

    /* initial misfit */
    cpxeiko_set(tr,ti);

/*
    cpxeiko_forw(false,ti,w);
    
    for (it=0; it < nt; it++) {
	if (known[it]) {
	    wi[it] = fix[it];
	} else {
	    wi[it] = w[it];
	    if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
	}
	
	wr[it] = s[it]+wi[it];
    }
    
    fastmarchcpx(tr,tr0,known,wr);
    fastmarchcpx(ti,ti0,known,wi);

    cpxeiko_set(tr,ti);
*/
    
    cpxeiko_forw(false,ti,w);

/*
    for (it=0; it < nt; it++) 
	w[it] = wi[it];
*/  
    for (it=0; it < nt; it++) 
	x0[it] = fix[it]-w[it];
    
    cpxeiko_forw(false,tr,rhsr);
    cpxeiko_forw(true, ti,rhsi);
    
    for (it=0; it < nt; it++) {
	rhs[it] = wt[it]*(-rhsr[it]-rhsi[it]);
    }
    
    if (NULL != witer) sf_floatwrite(w,nt,witer);
    if (NULL != rhsiter) sf_floatwrite(rhs,nt,rhsiter);
    if (NULL != x0iter) sf_floatwrite(x0,nt,x0iter);
    if (NULL != titer) {
	for (it=0; it < nt; it++)
	    t[it] = sf_cmplx(tr[it],ti[it]);
	
	sf_complexwrite(t,nt,titer);
    }

    /* calculate L2 data-misfit */
    rhsnorm0 = cblas_snrm2(nt,rhs,1);
    rhsnorm1 = rhsnorm0;
    rhsnorm = rhsnorm0;
    
    sf_warning("Iteration 0 of %d:\t misfit=%g, w=%g.",niter,rhsnorm,cblas_snrm2(nt,w,1));

    if (NULL != upiter) {
	cpxeiko_sten(pdir);
	sf_complexwrite(pdir[0],nt*dim,upiter);
    }

    if (NULL != operiter) {
	count = 0;

	for (i=0; i < nt; i++) {
	    for (it=0; it < nt; it++) {
		if (it == i)
		    unit[it] = 1.;
		else
		    unit[it] = 0.;
	    }

	    cpxeiko_operator(false,false,nt,nt,unit,oper[0][i]);
	    cpxeiko_operator(true,false,nt,nt,oper[1][i],unit);
	    
	    count++;
	}

	sf_floatwrite(oper[0][0],nt*nt*2,operiter);
    }

    for (iter=0; iter < niter; iter++) {

	/* clean-up */
	for (it=0; it < nt; it++) {
	    dw[it] = 0.;
	}
	
	/* solve dw */
	sf_solver(cpxeiko_operator,sf_cgstep,nt,nt,dw,rhs,cgiter,"known",known,"x0",x0,"wt",wt,"verb",verb,"end");
	sf_cgstep_close();
	
	if (NULL != dwiter) sf_floatwrite(dw,nt,dwiter);

	if (NULL != liniter) {
	    cpxeiko_operator(false,false,nt,nt,dw,rhsr);
	    sf_floatwrite(rhsr,nt,liniter);
	}

	gama = 1.;
	tol = 1.e-8;
	
        /* trying to avoid the points where w is very close to zero */
	for (it=0; it < nt; it++) {
	    if (dw[it] < 0. && w[it] > tol) {
		ratio = -w[it]/dw[it];
		gama = (gama<ratio)?gama:ratio;
	    }
	}

	/* linesearch */
	for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {

	    /* update real and imaginary slowness */
	    for (it=0; it < nt; it++) {
		if (known[it]) {
		    wi[it] = w[it]+dw[it];
		} else {
		    wi[it] = w[it]+gama*step*dw[it];
		    
		    /* clip negative slowness squared */
		    if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
		}
		
		wr[it] = s[it]+wi[it];
	    }
	    
	    /* forward fast-marching for stencil time */
	    fastmarchcpx(tr,tr0,known,wr);
	    fastmarchcpx(ti,ti0,known,wi);
	    
	    cpxeiko_set(tr,ti);
	    
	    cpxeiko_forw(false,tr,rhsr);
	    cpxeiko_forw(true, ti,rhsi);	    
	    for (it=0; it < nt; it++) {
		rhs[it] = wt[it]*(-rhsr[it]-rhsi[it]);
	    }
	    
	    rhsnorm = cblas_snrm2(nt,rhs,1);
	    
	    /* break */
	    if (rhsnorm < rhsnorm1) break;
	}	

/*
	if (istep == nstep) {
	    for (it=0; it < nt; it++) {
		if (known[it]) {
		    wi[it] = w[it]+dw[it];
		} else {
		    wi[it] = w[it];
		    if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
		}
		
		wr[it] = s[it]+wi[it];
	    }
	}

	fastmarchcpx(tr,tr0,known,wr);
	fastmarchcpx(ti,ti0,known,wi);
	
	cpxeiko_set(tr,ti);

	cpxeiko_forw(false,tr,rhsr);
	cpxeiko_forw(true, ti,rhsi);	    
	for (it=0; it < nt; it++) {
	    rhs[it] = wt[it]*(-rhsr[it]-rhsi[it]);
	}

	rhsnorm = cblas_snrm2(nt,rhs,1);
*/
	
	rhsnorm1 = rhsnorm;
/*	
	cpxeiko_forw(false,ti,w);
*/
	for (it=0; it < nt; it++) 
	    w[it] = wi[it];

	for (it=0; it < nt; it++) 
	    x0[it] = fix[it]-w[it];
  	
	if (NULL != witer) sf_floatwrite(w,nt,witer);
	if (NULL != rhsiter) sf_floatwrite(rhs,nt,rhsiter);
	if (NULL != x0iter) sf_floatwrite(x0,nt,x0iter);
	if (NULL != titer) {
	    for (it=0; it < nt; it++)
		t[it] = sf_cmplx(tr[it],ti[it]);
	    
	    sf_complexwrite(t,nt,titer);
	}
	if (NULL != wtiter) {
	    for (it=0; it < nt; it++)
		t[it] = sf_cmplx(wr[it],wi[it]);
	    
	    sf_complexwrite(t,nt,wtiter);
	}
	
	sf_warning("Iteration %d of %d:\t misfit=%g, w=%g, dw=%g (gama %g, line-search %d).",
		   iter+1,niter,rhsnorm,cblas_snrm2(nt,w,1),cblas_snrm2(nt,dw,1),gama,istep);
    }
    
    for (it=0; it < nt; it++)
	t[it] = sf_cmplx(tr[it],ti[it]);

    sf_complexwrite(t,nt,out);
  
    exit(0);
}
