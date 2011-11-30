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
    bool *knownr, *knowni, velocity, verb, *tknown, recom, pvar, reg, term, wupg, smooth;
    int dim, i, n[SF_MAX_DIM], *m, *tm, it, nt, iter, niter, cgiter, istep, nstep, rect[SF_MAX_DIM], irep, repeat;
    float d[SF_MAX_DIM], o[SF_MAX_DIM], *s, *tr, *ti, *tr0, *ti0, *gammat, *scale, *cost;
    float *w, *wr, *wi, *dw, *dws, *rhs, *rhsr, *rhsi, *x0, *fix, *wt, ***oper, ***mat1, ***mat2, *unit, *w0, *wir, *wri;
    float gama, ratio, tol, eps, namda, alpha;
    sf_complex *t, *t0, **pdir;
    float rhsnorm, rhsnorm0, rhsnorm1, step;
    char key[6], *prec, *bound, *symm;
    sf_file in, out, vel, maskr, maski, ref, wght, cray;
    sf_file witer, dwiter, dwsiter, rhsiter, upiter, x0iter, titer, wtiter, liniter, operiter, matriter, matiiter, gamiter, preciter;
    
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
    
    if (!sf_getfloat("tol",&tol)) tol=1.e-8;
    /* thresholding for gradient scaling */

    if (NULL == (symm = sf_getstring("symm"))) symm="both";
    /* right-hand side evaluation L_R*I or L_I*R (default both) */

    if(!sf_getbool("wupg",&wupg)) wupg=true;
    /* compute w for angle preconditioning */

    if (!sf_getbool("term",&term)) term=false;
    /* early termination if line-search failure */

    if (!sf_getbool("smooth",&smooth)) smooth=false;
    /* smooth update after conjugate-gradient */

    if (!sf_getint("repeat",&repeat)) repeat=1;
    /* number of smoothings */

    if (NULL == (prec = sf_getstring("prec"))) prec="angle";
    /* rhs preconditioning (default angle) */

    if (NULL == (bound = sf_getstring("bound"))) bound="add";
    /* avoid overshoot when update (default add) */

    if (!sf_getbool("reg",&reg)) reg=false;
    /* regularization (Ticknov) */

    if (!sf_getfloat("eps",&eps)) eps=1.e-2;
    /* stable division of preconditioner */
    
    if (!sf_getfloat("namda",&namda)) namda=0.1;
    /* regularization parameter (Ticknov) */

    if (!sf_getfloat("alpha",&alpha)) alpha=1.;
    /* exponential scaling of preconditioning */

    if (!sf_getbool("pvar",&pvar)) pvar=true;
    /* allow preconditioning to change over iterations */

    if (!sf_getbool("recom",&recom)) recom=true;
    /* recompute initial R according to w estimated from I */
    
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

    /* output dws at each iteration */
    if (NULL != sf_getstring("dwsiter")) {
	dwsiter = sf_output("dwsiter");
	sf_settype(dwsiter,SF_FLOAT);
	sf_putint(dwsiter,"n3",n[2]);
	sf_putint(dwsiter,"n4",niter);
    } else {
	dwsiter = NULL;
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
    } else {
	operiter = NULL;
	oper = NULL;
    }
    
    /* output operator at each iteration */
    /* NOTE: only first iteration */
    if (NULL != sf_getstring("matriter")) {
	matriter = sf_output("matriter");
	sf_settype(matriter,SF_FLOAT);
	sf_putint(matriter,"n1",nt);
	sf_putint(matriter,"n2",nt);
	sf_putint(matriter,"n3",4);
	sf_putint(matriter,"n4",1);
	mat1 = sf_floatalloc3(nt,nt,4);
    } else {
	matriter = NULL;
	mat1 = NULL;
    }
    
    /* output operator at each iteration */
    /* NOTE: only first iteration */
    if (NULL != sf_getstring("matiiter")) {
	matiiter = sf_output("matiiter");
	sf_settype(matiiter,SF_FLOAT);
	sf_putint(matiiter,"n1",nt);
	sf_putint(matiiter,"n2",nt);
	sf_putint(matiiter,"n3",4);
	sf_putint(matiiter,"n4",1);
	mat2 = sf_floatalloc3(nt,nt,4);
    } else {
	matiiter = NULL;
	mat2 = NULL;
    }
    
    /* unit vector */
    if (operiter != NULL || matriter != NULL || matiiter != NULL) {
	unit = sf_floatalloc(nt);
    } else {
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
    
    /* output linesearch at 1st iteration */
    if (NULL != sf_getstring("liniter")) {
	liniter = sf_output("liniter");
	sf_settype(liniter,SF_FLOAT);
	sf_putint(liniter,"n3",n[2]);
	sf_putint(liniter,"n4",nstep);
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
    
    /* output gamma at each iteration */
    if (NULL != sf_getstring("gamiter")) {
	gamiter = sf_output("gamiter");
	sf_settype(gamiter,SF_FLOAT);
	sf_putint(gamiter,"n3",n[2]);
	sf_putint(gamiter,"n4",niter);
	gammat = sf_floatalloc(nt);
    } else {
	gamiter = NULL;
	gammat = NULL;
    }
    
    /* output preconditioning at each iteration */
    if (NULL != sf_getstring("preciter")) {
	preciter = sf_output("preciter");
	sf_settype(preciter,SF_FLOAT);
	sf_putint(preciter,"n3",n[2]);
	sf_putint(preciter,"n4",niter+1);
    } else {
	preciter = NULL;
    }
    
    /* read boundary condition */
    m = sf_intalloc(nt);
    knownr = sf_boolalloc(nt);
    knowni = sf_boolalloc(nt);

    if (NULL != sf_getstring("maskr")) {
	maskr = sf_input("maskr");
	sf_intread(m,nt,maskr);
	sf_fileclose(maskr);
	
	for (it=0; it < nt; it++) {
	    if (m[it] != 0)
		knownr[it] = true;
	    else
		knownr[it] = false;
	}
    } else {
	for (it=0; it < nt; it++) {
	    knownr[it] = false;
	}
    }
    
    if (NULL != sf_getstring("maski")) {
	maski = sf_input("maski");
	sf_intread(m,nt,maski);
	sf_fileclose(maski);
	
	for (it=0; it < nt; it++) {
	    if (m[it] != 0)
		knowni[it] = true;
	    else
		knowni[it] = false;
	}
    } else {
	for (it=0; it < nt; it++) {
	    knowni[it] = false;
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
    
    /* NOTE: w_exact */
    cpxeiko_ref(dim,n,d,ti0,fix);
    
    free(m);
    free(t0);
    m = NULL;
    t0 = NULL;
    
    /* read right-hand side weighting */
    wt = sf_floatalloc(nt);
    scale = sf_floatalloc(nt);
    cost = sf_floatalloc(nt);
    
    if (NULL != sf_getstring("wght")) {
	wght = sf_input("wght");
	sf_floatread(scale,nt,wght);
	sf_fileclose(wght);
    } else {
	for (it=0; it < nt; it++) {
	    scale[it] = 1.;
	}
    }
    
    /* allocate temporary memory */
    w    = sf_floatalloc(nt);
    rhsr = sf_floatalloc(nt);
    rhsi = sf_floatalloc(nt);
    wr   = sf_floatalloc(nt);
    wi   = sf_floatalloc(nt);
    x0   = sf_floatalloc(nt);
    dw   = sf_floatalloc(nt);
    dws  = sf_floatalloc(nt);
    rhs  = sf_floatalloc(nt);
    
    /* initialize 2D gradient operator */
    sf_igrad2_init(n[0],n[1]);
    
    /* initialize triangular smoother */
    for (i=0; i < dim; i++) {
	sprintf(key,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/
    }
    sf_trianglen_init(dim,rect,n);

    /* initialize fastmarchcpx */
    /* NOTE: default accuracy 1st order */
    fastmarchcpx_init(n,o,d,1);
    
    /* initialize cpxeiko */
    cpxeiko_init(dim,n,nt,d);
    
    /* initial misfit */
    cpxeiko_set(tr,ti);
    
    /* w Gaussian beam */
    w0 = sf_floatalloc(nt);
    cpxeiko_forw(false,ti,w0);

/* NOTE: the following lines recompute initial R and I according to w from I. */
    if (recom) {
	tknown = sf_boolalloc(nt);

	/* NOTE: w_gaussian */
	cpxeiko_forw(false,ti,w);
	
	if (NULL != sf_getstring("cray")) {
	    cray = sf_input("cray");
	    tm = sf_intalloc(nt);
	    sf_intread(tm,nt,cray);
	    sf_fileclose(cray);

	    for (it=0; it < nt; it++) {
		if (tm[it] == 1) {
		    tknown[it] = true;
		} else {
		    tknown[it] = false;
		}

		wi[it] = w[it];
		wr[it] = s[it]+wi[it];
	    }
		
	    free(tm);
	} else {
	    for (it=0; it < nt; it++) {
		if (ti[it] < 1.e-8) {
		    tknown[it] = true;
		} else {
		    tknown[it] = false;
		}
		
		wi[it] = w[it];
		wr[it] = s[it]+wi[it];
	    }
	}
	
	fastmarchcpx(tr,tr0,tknown,wr);
	fastmarchcpx(ti,ti0,tknown,wi);
	
	cpxeiko_set(tr,ti);
    } else {
	tknown = NULL;

	cpxeiko_forw(false,ti,wi);
    }
    
    for (it=0; it < nt; it++) 
	w[it] = wi[it];
    
    /* NOTE: replace the central ray tube of fix with w Gaussian beam (2D only) */
    for (i=0; i < n[1]; i++) {
	for (it=1; it < n[0]; it++) {
	    if (knownr[i*n[0]+it]) {
		fix[i*n[0]+it] = w[i*n[0]+it];
		tr0[i*n[0]+it] = tr[i*n[0]+it];
		ti0[i*n[0]+it] = ti[i*n[0]+it];
	    }
	}
    }

    /* right-hand side preconditioning */
    /* NOTE: for ray-traced beam the central ray may not be charaterized by w0[it]==0. 
       use eps != 0. and scale == 0. to avoid 0./0. = NAN issue. */
    if (wupg) {
	wir = sf_floatalloc(nt);
	wri = sf_floatalloc(nt);
    } else {
	wir = NULL;
	wri = NULL;
    }

    switch (prec[0]) {
	case 'n': /* none */

	    for (it=0; it < nt; it++) {
		wt[it] = scale[it];
	    }

	    break;

	case 'a': /* angle */

	    if (wupg) {
		cpxeiko_wupg(false,tr,wir);
		cpxeiko_wupg(true,ti,wri);
	    }

	    for (it=0; it < nt; it++) {
		if (pvar) {
		    if (wupg) {
			/*
			  if (wi[it] > 0.)
			*/
			if (w0[it] > 0.)
			    wt[it] = scale[it]/(sqrtf(powf(wi[it],alpha)*wir[it])+eps);
			else
			    wt[it] = 0.;
		    } else {
			/*
			  if (wi[it] > 0.)
			*/
			if (w0[it] > 0.)
			    wt[it] = scale[it]/(sqrtf(powf(wi[it],alpha)*(s[it]+wi[it]))+eps);
			else
			    wt[it] = 0.;
		    }
		} else {
		    if (w0[it] > 0.)
			wt[it] = scale[it]/(sqrtf(powf(w0[it],alpha)*(s[it]+w0[it]))+eps);
		    else
			wt[it] = 0.;
		}
	    }
	    
	    break;
	    
	case 'e': /* exponential */
	    
	    for (it=0; it < nt; it++) {
		if (pvar) {
		    /*
		      if (wi[it] > 0.)
		    */
		    if (w0[it] > 0.)
			wt[it] = scale[it]/(expf(wi[it])+eps);
		    else
			wt[it] = 0.;
		} else {
		    if (w0[it] > 0.)
			wt[it] = scale[it]/(expf(w0[it])+eps);
		    else
			wt[it] = 0.;
		}
	    }
	    
	    break;
    }
    
    for (it=0; it < nt; it++) 
	x0[it] = fix[it]-w[it];
    
    cpxeiko_forw(false,tr,rhsr);
    cpxeiko_forw(true, ti,rhsi);
    
    for (it=0; it < nt; it++) {
	switch (symm[0]) {
	    case 'b': /* both L_R*I and L_I*R */
		rhs[it] = -rhsr[it]-rhsi[it];
		break;

	    case 'r': /* L_R*I */
		rhs[it] = -rhsi[it]-rhsi[it];
		break;

	    case 'i': /* L_I*R */
		rhs[it] = -rhsr[it]-rhsr[it];
		break;
	}
	
	cost[it] = wt[it]*rhs[it];
    }
    
    if (NULL != witer) sf_floatwrite(w,nt,witer);
    if (NULL != rhsiter) sf_floatwrite(rhs,nt,rhsiter);
    if (NULL != x0iter) sf_floatwrite(x0,nt,x0iter);
    if (NULL != titer) {
	for (it=0; it < nt; it++)
	    t[it] = sf_cmplx(tr[it],ti[it]);
	
	sf_complexwrite(t,nt,titer);
    }
    if (NULL != preciter) sf_floatwrite(wt,nt,preciter);

    /* calculate L2 data-misfit */
    rhsnorm0 = cblas_snrm2(nt,cost,1);
    rhsnorm1 = rhsnorm0;
    rhsnorm = rhsnorm0;
    
    sf_warning("Iteration 0 of %d:\t misfit=%g, w=%g.",niter,rhsnorm,cblas_snrm2(nt,w,1));

    if (NULL != upiter) {
	cpxeiko_sten(pdir);
	sf_complexwrite(pdir[0],nt*dim,upiter);
    }

    /* output operators in matrix form */
    if (NULL != operiter || NULL != matriter || NULL != matiiter ) {
	for (i=0; i < nt; i++) {
	    for (it=0; it < nt; it++) {
		if (it == i)
		    unit[it] = 1.;
		else
		    unit[it] = 0.;
	    }

	    if (NULL != operiter) {
		cpxeiko_operator(false,false,nt,nt,unit,oper[0][i]);
		cpxeiko_operator(true,false,nt,nt,oper[1][i],unit);
	    }	    

	    if (NULL != matriter) {
		cpxeiko_mat(false,0,nt,nt,unit,mat1[0][i]);
		cpxeiko_mat(false,1,nt,nt,unit,mat1[1][i]);
		cpxeiko_mat(false,2,nt,nt,unit,mat1[2][i]);
		cpxeiko_mat(false,3,nt,nt,unit,mat1[3][i]);
	    }
	    if (NULL != matiiter) {
		cpxeiko_mat(true,0,nt,nt,unit,mat2[0][i]);
		cpxeiko_mat(true,1,nt,nt,unit,mat2[1][i]);
		cpxeiko_mat(true,2,nt,nt,unit,mat2[2][i]);
		cpxeiko_mat(true,3,nt,nt,unit,mat2[3][i]);
	    }
	}
	
	if (NULL != operiter) {
	    sf_floatwrite(oper[0][0],nt*nt*2,operiter);
	}
	if (NULL != matriter) {
	    sf_floatwrite(mat1[0][0],nt*nt*4,matriter);
	}
	if (NULL != matiiter) {
	    sf_floatwrite(mat2[0][0],nt*nt*4,matiiter);
	}
    }

    for (iter=0; iter < niter; iter++) {

	/* clean-up */
	for (it=0; it < nt; it++) {
	    dw[it] = 0.;
	}

	/* solve dw */
	if (reg)
	    sf_solver_reg(cpxeiko_operator,sf_cgstep,sf_igrad2_lop,2*nt,nt,nt,dw,rhs,cgiter,namda,"known",knowni,"x0",x0,"wt",wt,"verb",verb,"end");
	else
	    sf_solver(cpxeiko_operator,sf_cgstep,nt,nt,dw,rhs,cgiter,"known",knownr,"x0",x0,"wt",wt,"verb",verb,"end");
	
	sf_cgstep_close();
	
	if (NULL != dwiter) sf_floatwrite(dw,nt,dwiter);
	
	/* NOTE: smooth dw before line-search */
/*
	if (smooth) {
	    for (irep=0; irep < repeat; irep++) {
		sf_trianglen_lop(false,false,nt,nt,dw,dws);

		for (it=0; it < nt; it++) {
		    if (knownr[it])
			dw[it] = x0[it];
		    else
			dw[it] = dws[it];
		}
	    }
	}

	if (NULL != dwsiter) sf_floatwrite(dw,nt,dwsiter);
*/
	gama = 1.;
	
        /* trying to avoid the points where w is very close to zero */
	for (it=0; it < nt; it++) {
	    if (dw[it] < 0. && w[it] > tol) {
		ratio = -w[it]/dw[it];

		if (gama > ratio) gama = ratio;
	    }
	    
	    if (NULL != gammat) {
		if (w[it]+dw[it] < 0.)
		    gammat[it] = 1.;
		else
		    gammat[it] = 0.;
	    }
	}
	
	if (NULL != gamiter) sf_floatwrite(gammat,nt,gamiter);

	/* linesearch */
	for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {

	    /* update real and imaginary slowness */
	    for (it=0; it < nt; it++) {
		if (knowni[it]) {
		    wi[it] = w[it]+dw[it];
		} else {
		    switch (bound[0]) {
			case 'a': /* add */

			    wi[it] = w[it]+gama*step*dw[it];
			    
			    break;

			case 'e': /* exponential */

			    wi[it] = w[it]*expf(gama*step*dw[it]/w[it]);

			    break;

			case 'h': /* hybrid */

			    /* DEBUG!!! */
			    if (dw[it] < 0. && w[it] <= namda)
				wi[it] = w[it]*expf(gama*step*dw[it]/w[it]);
			    else
				wi[it] = w[it]+gama*step*dw[it];

			    break;
		    }
		    
		    /* NOTE: what to do with overshoot (lower bound) */
		    /*
		    if (wi[it] <= 0.) wi[it] = FLT_EPSILON;
		    */
		    /*
		    if (wi[it] <= 0.) wi[it] = w[it];
		    */
		    if (wi[it] <= 0.) wi[it] = 0.;

		    /* NOTE: what to do with overshoot (upper bound) */
		    /*
		      
		    */
		    }
		
		wr[it] = s[it]+wi[it];
	    }
	    
	    if (smooth) {
		for (irep=0; irep < repeat; irep++) {
		    sf_trianglen_lop(false,false,nt,nt,wi,dws);
		    for (it=0; it < nt; it++) {
			if (knowni[it])
			    wi[it] = w[it]+dw[it];
			else
			    wi[it] = dws[it];
		    }
		    
		    sf_trianglen_lop(false,false,nt,nt,wr,dws);
		    for (it=0; it < nt; it++) {
			if (knowni[it])
			    wr[it] = s[it]+wi[it];
			else
			    wr[it] = dws[it];
		    }
		}
	    }

	    /* forward fast-marching for stencil time */

	    fastmarchcpx(tr,tr0,knownr,wr);
	    /* NOTE: only supply I boundary condition on central ray
	    fastmarchcpx(ti,ti0,knowni,wi);
	    */
	    fastmarchcpx(ti,ti0,knowni,wi);

	    /* right-hand side preconditioning */
	    /* NOTE: for ray-traced beam the central ray may not be charaterized by w0[it]==0. 
	       use eps != 0. and scale == 0. to avoid 0./0. = NAN issue. */
	    switch (prec[0]) {
		case 'n': /* none */
		    
		    for (it=0; it < nt; it++) {
			wt[it] = scale[it];
		    }
		    
		    break;
		    
		case 'a': /* angle */
		    
		    if (wupg) {
			cpxeiko_wupg(false,tr,wir);
			cpxeiko_wupg(true,ti,wri);
		    }
		    
		    for (it=0; it < nt; it++) {
			if (pvar) {
			    if (wupg) {
				/*
				  if (wi[it] > 0.)
				*/
				if (w0[it] > 0.)
				    wt[it] = scale[it]/(sqrtf(powf(wi[it],alpha)*wir[it])+eps);
				else
				    wt[it] = 0.;
			    } else {
				/*
				  if (wi[it] > 0.)
				*/
				if (w0[it] > 0.)
				    wt[it] = scale[it]/(sqrtf(powf(wi[it],alpha)*(s[it]+wi[it]))+eps);
				else
				    wt[it] = 0.;
			    }
			} else {
			    if (w0[it] > 0.)
				wt[it] = scale[it]/(sqrtf(powf(w0[it],alpha)*(s[it]+w0[it]))+eps);
			    else
				wt[it] = 0.;
			}
		    }
		    
		    break;
		    
		case 'e': /* exponential */
		    
		    for (it=0; it < nt; it++) {
			if (pvar) {
			    /*
			      if (wi[it] > 0.)
			    */
			    if (w0[it] > 0.)
				wt[it] = scale[it]/(expf(wi[it])+eps);
			    else
				wt[it] = 0.;
			} else {
			    if (w0[it] > 0.)
				wt[it] = scale[it]/(expf(w0[it])+eps);
			    else
				wt[it] = 0.;
			}
		    }

		    break;
	    }
	    
	    cpxeiko_set(tr,ti);
	    
	    cpxeiko_forw(false,tr,rhsr);
	    cpxeiko_forw(true, ti,rhsi);	    
	    for (it=0; it < nt; it++) {
		switch (symm[0]) {
		    case 'b': /* both L_R*I and L_I*R */
			rhs[it] = -rhsr[it]-rhsi[it];
			break;
			
		    case 'r': /* L_R*I */
			rhs[it] = -rhsi[it]-rhsi[it];
			break;
			
		    case 'i': /* L_I*R */
			rhs[it] = -rhsr[it]-rhsr[it];
			break;
		}
		
		cost[it] = wt[it]*rhs[it];
	    }
	    
	    if (iter == 0 && NULL != liniter) {
		sf_floatwrite(cost,nt,liniter);
	    }
	    
	    rhsnorm = cblas_snrm2(nt,cost,1);
	    
	    /* break */
	    if (rhsnorm < rhsnorm1) break;
	}	

	/* early termination if line-search failed */
	if (term && (istep == nstep)) {
	    sf_warning("Line-search failure. Output dimensions need to be fixed before read.");

	    break;
	}

/* NOTE: the following lines supplies boundary condition but force 
   other w to be the same as previous iteration. */
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
/* NOTE: the following line is trying to recompute w from forward modeled I.
   This will introduce differences at least at boundaries. */
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
	if (NULL != preciter) sf_floatwrite(wt,nt,preciter);
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
