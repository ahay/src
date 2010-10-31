/* Linearized complex eikonal equation */
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

#include "cpxeiko.h"
#include "fastmarchcpx.h"

int main(int argc, char* argv[])
{
    bool *m, verb;
    int dim, i, n[SF_MAX_DIM], it, nt, iter, niter, cgiter, *temp, is, ns, *intm;
    float d[SF_MAX_DIM], o[SF_MAX_DIM], *tr, *ti, *tr0, *ti0, *w, *dw, *s, *rhs, *x0, *wi, *wr, *fix, *exact;
    float dwnorm=0.0, wnorm0=0.0, rhsnorm=0.0, rhsnorm0=1.0, rhsnorm2=2.0, step;
/*    
    double dot1[2], dot2[2]; 
*/
    char key[4];
    sf_file realtime, imagtime, realout, imagout, slow, ref, monitor, slowout;
    
    sf_init(argc,argv);
    realtime = sf_input("in");
    imagtime = sf_input("imagtime");
    slow = sf_input("slow");
    realout = sf_output("out");
    imagout = sf_output("imagout");

    /* read in model dimensions */
    dim = sf_filedims(realtime,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(realtime,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(realtime,key,o+i)) sf_error("No %s= in input",key);
	nt *= n[i];
    }

    if (dim < 3) {
	n[2] = 1; d[2] = d[1]; o[2] = o[1];
    }

    if(!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
    /* number of conjugate gradient iterations */

    if (!sf_getint("ns",&ns)) ns=10;
    /* number of iterations in step size search */

    /* optional output: slowness */
    if (NULL != sf_getstring("slowout")) {
	slowout = sf_output("slowout");
	sf_putint(slowout,"n3",n[2]);
	sf_putint(slowout,"n4",niter);
    } else {
	slowout = NULL;
    }

    /* optional output: L2 misfit */
    if (NULL != sf_getstring("monitor")) {
	monitor = sf_output("monitor");
	sf_putint(monitor,"n3",n[2]);
	sf_putint(monitor,"n4",niter);
    } else {
	monitor = NULL;
    }

    m     = sf_boolalloc(nt);
    w     = sf_floatalloc(nt);
    s     = sf_floatalloc(nt);
    tr    = sf_floatalloc(nt);
    ti    = sf_floatalloc(nt);
    x0    = sf_floatalloc(nt);
    dw    = sf_floatalloc(nt);
    rhs   = sf_floatalloc(nt);
    fix   = sf_floatalloc(nt);
    temp  = sf_intalloc(nt);
    exact = sf_floatalloc(nt);
    intm  = sf_intalloc(nt);
    wi    = sf_floatalloc(nt);
    wr    = sf_floatalloc(nt);
/*
    if (NULL != sf_getstring("mask"))
    {
	mask = sf_input("mask");
	if (SF_INT != sf_gettype(mask))
	    sf_error("need int mask input");
	sf_intread(intm,nt,mask);
	sf_fileclose(mask);
	
	for (it=0; it < nt; it++) {
	    if (intm[it] == 1)
		m[it] = true;
	    else
		m[it] = false;
	}
    }
*/
    /* fix later for 3-D */

    ti0 = sf_floatalloc(n[1]);
    tr0 = sf_floatalloc(n[1]);

    /* input exact solution */   
    ref = sf_input("realref");
    sf_floatread(exact,nt,ref);
    sf_fileclose(ref);

    for (it=0; it < n[1]; it++) {
	tr0[it] = exact[it*n[0]];
    }

    ref = sf_input("imagref");
    sf_floatread(exact,nt,ref);
    sf_fileclose(ref);

    for (it=0; it < n[1]; it++) {
	ti0[it] = exact[it*n[0]];
    }

    /* compute exact w */
    cpxeiko_ref(dim,n,d,exact,fix);

    sf_floatread(tr,nt,realtime);
    sf_floatread(ti,nt,imagtime);
    sf_floatread(s,nt,slow);

    fastmarchcpx_init(n,d);
    cpxeiko_init(dim,n,nt,d);

    for (iter=0; iter < niter; iter++) {
	if (iter == 0) {
	    /* initiate operator */
	    cpxeiko_set(tr,ti);
	    
	    /* current pseudo slowness */
	    cpxeiko_forw(false,ti,w);

	    for (it=0; it < nt; it++) {
		m[it] = false;
/*		x0[it]  = 0.0; */
	    }

	    /* first two layers */
	    for (it=0; it < n[1]; it++) {
/*		
                w[it*n[0]]   = fix[it*n[0]];
		w[it*n[0]+1] = fix[it*n[0]+1];
*/		
		m[it*n[0]]   = true;
	    }

	    wnorm0 = cblas_snrm2(nt, w, 1);

	    /*

	    for (it=0; it < nt; it++) {
		wr[it] = (1/s[it])*(1/s[it])+w[it];
	    }

	    fastmarchcpx(tr,tr0,wr);

	    cpxeiko_set(tr,ti);

	    */

	    cpxeiko_forw(false,tr,wr);
	    cpxeiko_forw(true, ti,wi);

	    /* right-hand side */	
	    for (it=0; it < nt; it++) {
		rhs[it] = -wr[it]-wi[it]; 
	    }

	    rhsnorm0 = cblas_snrm2(nt, rhs, 1);
	    rhsnorm  = 1.0;
	    rhsnorm2 = 2.0;
	}

	for (it=0; it < nt; it++) {
	    x0[it] = fix[it]-w[it];
	}
/*
	if (NULL != slowout) sf_floatwrite(w,nt,slowout);
	if (NULL != monitor) sf_floatwrite(rhs,nt,monitor);
*/

	if (NULL != slowout) sf_floatwrite(w,nt,slowout);

	sf_solver(cpxeiko_operator,sf_cgstep,nt,nt,dw,rhs,cgiter,"known",m,"x0",x0,"verb",false,"end");
	sf_cgstep_close();

	if (NULL != monitor) sf_floatwrite(dw,nt,monitor);	

/*
  init_genrand((unsigned long) time(NULL));
  
  sf_dot_test(cpxeiko_operator,nt,nt,dot1,dot2);
  sf_warning("compare %g and %g",dot1[0],dot1[1]);
  sf_warning("compare %g and %g",dot2[0],dot2[1]);
*/

	/* loop over step size */
	for (is=0, step=1.0; is < ns; is++, step *= 0.5) {	
	    for (it=0; it < nt; it++) {
		wi[it] = w[it] + step*dw[it];

		if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
 	    
		wr[it] = (1/s[it])*(1/s[it])+wi[it];
	    }
	    
	    fastmarchcpx(ti,ti0,wi);
	    fastmarchcpx(tr,tr0,wr);

	    cpxeiko_set(tr,ti);

	    cpxeiko_forw(false,tr,wr);
	    cpxeiko_forw(true, ti,wi);

	    for (it=0; it < nt; it++) {
		rhs[it] = -wr[it]-wi[it];
	    }

	    rhsnorm2 = cblas_snrm2(nt, rhs, 1)/rhsnorm0;
	    if (rhsnorm2 < rhsnorm) break;
	}	

	rhsnorm = rhsnorm2;
	for (it=0; it < nt; it++) {
	    w[it] += step*dw[it];
	    
	    if (w[it] <= 0.) w[it] = FLT_EPSILON;
	}

	if (verb) {
	    dwnorm = cblas_snrm2(nt, dw, 1)/wnorm0;
	    sf_warning("Iteration: %d of %d (%g %g %d).",iter+1,niter,dwnorm,rhsnorm,is);
	}
    }

    sf_floatwrite(tr,nt,realout);
    sf_floatwrite(ti,nt,imagout);
  
    exit(0);
}
