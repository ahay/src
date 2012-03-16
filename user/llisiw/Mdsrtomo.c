/* Prestack first-arrival traveltime tomography (DSR) */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "dsreiko.h"
#include "dsrtomo.h"

int main(int argc, char* argv[])
{
    bool velocity, verb, adj, shape;
    int dimw, dimt, i, j, k, n[SF_MAX_DIM], rect[SF_MAX_DIM], iw, nw, it, nt;
    int iter, niter, cgiter, count;
    int *f, *m0, offset;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *dt, *dw, *dv=NULL, *t, *w, *t0, *w1, *p=NULL;
    float eps, tol, tau1, tau2, angle, *th, rhsnorm, rhsnorm0, rhsnorm1, rate, gama, *den=NULL;
    char key[6], *what;
    sf_file in, out, reco, grad, flag, mask, debug;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="tomo";
    /* what to compute (default tomography) */

    switch (what[0]) {
	case 'l': /* linear operator */

	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */

	    /* read record file */
	    if (NULL == sf_getstring("reco"))
		sf_error("Need time reco=");
	    reco = sf_input("reco");

	    /* read velocity file */
	    if (NULL == sf_getstring("grad"))
		sf_error("Need velocity grad=");
	    grad = sf_input("grad");
	    
	    /* read flag file */
	    if (NULL == sf_getstring("flag"))
		flag = NULL;
	    else
		flag = sf_input("flag");

	    /* read receiver file */
	    if (NULL == sf_getstring("mask"))
		mask = NULL;
	    else
		mask = sf_input("mask");

	    /* output debug file */
	    if (NULL == sf_getstring("debug"))
		debug = NULL;
	    else
		debug = sf_output("debug");	    

	    /* read dimension */
	    dimt = sf_filedims(reco,n);

	    nt = 1;
	    for (i=0; i < dimt; i++) {
		sprintf(key,"d%d",i+1);
		if (!sf_histfloat(reco,key,d+i)) sf_error("No %s= in input.",key);
		sprintf(key,"o%d",i+1);
		if (!sf_histfloat(reco,key,o+i)) o[i]=0.;
		nt *= n[i];
	    }

	    nw = nt/n[2];
	    
	    /* allocate memory */
	    w  = sf_floatalloc(nw);
	    dw = sf_floatalloc(nw);
	    t  = sf_floatalloc(nt);
	    dt = sf_floatalloc(nt);

	    if (flag != NULL)
		f  = sf_intalloc(nt);
	    else
		f = NULL;

	    m0  = sf_intalloc(nt);

	    /* read file */
	    sf_floatread(w,nw,grad);
	    sf_floatread(t,nt,reco);

	    if (flag != NULL) sf_intread(f,nt,flag);
	    if (mask != NULL) {
		sf_intread(m0,nt,mask);
	    } else {
		if (!sf_getint("offset",&offset)) offset=n[1];
		/* offset */

		for (k=0; k < n[2]; k++) {
		    for (j=0; j < n[1]; j++) {
			if (abs(j-k) <= offset)
			    m0[j*n[0]+k*n[0]*n[1]] = 1;
			else
			    m0[j*n[0]+k*n[0]*n[1]] = 0;

			for (i=1; i < n[0]; i++){
			    m0[i+j*n[0]+k*n[0]*n[1]] = 0;
			}
		    }
		}
	    }

	    if (debug != NULL) {
		den = sf_floatalloc(nt);

		if (!adj) {
		    sf_putint(debug,"n3",n[2]);
		    sf_putfloat(debug,"d3",d[2]);
		    sf_putfloat(debug,"o3",o[2]);
		}
	    }

	    if (!sf_getbool("velocity",&velocity)) velocity=true;
	    /* if y, the input is velocity; n, slowness squared */
	    
	    /* convert to slowness squared */
	    if (velocity) {
		for (iw=0; iw < nw; iw++)
		    w[iw] = 1./w[iw]*1./w[iw];

		dv = sf_floatalloc(nw);
	    }

	    if (adj) {
		sf_floatread(dt,nt,in);

		sf_putint(out,"n3",1);
	    } else {
		sf_floatread(dw,nw,in);

		sf_putint(out,"n3",n[2]);
		sf_putfloat(out,"d3",d[2]);
		sf_putfloat(out,"o3",o[2]);
	    }
	    
	    /* initialize operator */
	    dsrtomo_init(dimt,n,d);

	    /* set operator */
	    dsrtomo_set(t,w,f,m0);	    

	    if (adj) {
		dsrtomo_oper(true,false,nw,nt,dw,dt);
	    } else {
		dsrtomo_oper(false,false,nw,nt,dw,dt);	    
	    }
	    
	    if (debug != NULL) {
		dsrtomo_debug(den);

		sf_floatwrite(den,nt,debug);
	    }

	    if (adj) {
		/*
		if (velocity) {
		    for (iw=0; iw < nw; iw++)
			dv[iw] = -dw[iw]/(2.*sqrtf(w[iw])*(w[iw]+dw[iw]/2.));

		    sf_floatwrite(dv,nw,out);
		} else {
		    sf_floatwrite(dw,nw,out);
		}
		*/

		sf_floatwrite(dw,nw,out);
	    } else {
		sf_floatwrite(dt,nt,out);
	    }

	    break;
	    
	case 't': /* tomography */

	    /* read dimension */
	    dimw = sf_filedims(in,n);
	    
	    nw = 1;
	    for (i=0; i < dimw; i++) {
		sprintf(key,"d%d",i+1);
		if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
		sprintf(key,"o%d",i+1);
		if (!sf_histfloat(in,key,o+i)) o[i]=0.;
		nw *= n[i];
	    }
	    
	    if (dimw > 2) sf_error("Only works for 2D now.");
	    
	    n[2] = n[1]; d[2] = d[1]; o[2] = o[1]; 
	    dimt = 3; nt = nw*n[2];
	    
	    /* read initial velocity */
	    w = sf_floatalloc(nw);
	    sf_floatread(w,nw,in);
	    
	    if (!sf_getbool("velocity",&velocity)) velocity=true;
	    /* if y, the input is velocity; n, slowness squared */
	    
	    if (!sf_getbool("shape",&shape)) shape=false;
	    /* shaping regularization (default no) */

	    /* convert to slowness squared */
	    if (velocity) {
		for (iw=0; iw < nw; iw++)
		    w[iw] = 1./w[iw]*1./w[iw];

		dv = sf_floatalloc(nw);
	    }
	    
	    /* read record file */
	    if (NULL == sf_getstring("reco"))
		sf_error("Need record reco=");
	    reco = sf_input("reco");
	    
	    t0 = sf_floatalloc(nt);
	    sf_floatread(t0,nt,reco);
	    sf_fileclose(reco);
	    
	    /* read receiver file */
	    m0 = sf_intalloc(nt);

	    if (NULL == sf_getstring("mask")) {
		mask = NULL;
		
		if (!sf_getint("offset",&offset)) offset=n[1];
		/* offset */

		for (k=0; k < n[2]; k++) {
		    for (j=0; j < n[1]; j++) {
			if (abs(j-k) <= offset)
			    m0[j*n[0]+k*n[0]*n[1]] = 1;
			else
			    m0[j*n[0]+k*n[0]*n[1]] = 0;
			
			for (i=1; i < n[0]; i++){
			    m0[i+j*n[0]+k*n[0]*n[1]] = 0;
			}
		    }
		}
	    } else {
		mask = sf_input("mask");
		
		sf_intread(m0,nt,mask);
		sf_fileclose(mask);
	    }	    

	    if (!sf_getbool("verb",&verb)) verb=false;
	    /* verbosity flag */
	    
	    if (!sf_getint("niter",&niter)) niter=5;
	    /* number of inversion iterations */
	    
	    if (!sf_getint("cgiter",&cgiter)) cgiter=25;
	    /* number of conjugate-gradient iterations */
	    
	    if (!sf_getfloat("tau1",&tau1)) tau1=1.e-3;
	    /* tau1 */

	    if (!sf_getfloat("tau2",&tau2)) tau2=1.;
	    /* tau2 */

	    if (!sf_getfloat("angle",&angle)) angle=0.1;
	    /* angle */

	    /* output gradient at each iteration */
	    if (NULL != sf_getstring("grad")) {
		grad = sf_output("grad");
		sf_putint(grad,"n3",niter);
	    } else {
		grad = NULL;
	    }
	    
	    if (!sf_getfloat("eps",&eps)) eps=0.;
	    /* regularization parameter */

	    if (shape) {
		if (!sf_getfloat("tol",&tol)) tol=1.e-6;
		/* tolerance for shaping regularization */

		for (i=0; i < dimw; i++) {
		    sprintf(key,"rect%d",i+1);
		    if (!sf_getint(key,rect+i)) rect[i]=1;
		    /*( rect#=(1,1,...) smoothing radius on #-th axis )*/
		}
		
		/* triangle smoothing operator */
		sf_trianglen_init(dimw,rect,n);
		sf_repeat_init(nw,1,sf_trianglen_lop);
		
		sf_conjgrad_init(nw,nw,nt,nt,eps,tol,verb,false);
		p = sf_floatalloc(nw);
	    } else {
		/* initialize 2D gradient operator */
		sf_igrad2_init(n[0],n[1]);
	    }

	    /* allocate temporary array */
	    t  = sf_floatalloc(nt);
	    dw = sf_floatalloc(nw);
	    dt = sf_floatalloc(nt);
	    w1 = sf_floatalloc(nw);
	    f  = sf_intalloc(nt);
	    th = sf_floatalloc(nt);

	    /* initialize eikonal */
	    dsreiko_init(n,o,d,tau1,tau2,angle);
	    
	    /* initialize operator */
	    dsrtomo_init(dimt,n,d);
	    
	    /* initial misfit */
	    dsreiko_fastmarch(t,w,f,th);
	    dsreiko_mirror(t);
	    
	    /* calculate L2 data-misfit */
	    for (it=0; it < nt; it++) {
		if (m0 == NULL || m0[it] == 1)
		    dt[it] = t0[it]-t[it];
		else
		    dt[it] = 0.;
	    }
	    
	    rhsnorm0 = cblas_snrm2(nt,dt,1);
	    rhsnorm = rhsnorm0;
	    rhsnorm1 = rhsnorm;
	    rate = rhsnorm1/rhsnorm0;
	    
	    sf_warning("L2 misfit after iteration 0 of %d: %g",niter,rate);
	    
	    /* iterations over inversion */
	    for (iter=0; iter < niter; iter++) {
		
		/* clean-up */
		for (iw=0; iw < nw; iw++)
		    dw[iw] = 0.;
		
		/* set operator */
		dsrtomo_set(t,w,f,m0);
		
		/* solve dw */
		if (shape) {
		    sf_conjgrad(NULL,dsrtomo_oper,sf_repeat_lop,p,dw,dt,cgiter);
		} else {
		    /*
		    sf_solver(dsrtomo_oper,sf_cgstep,nw,nt,dw,dt,cgiter,"verb",verb,"end");
		    */
		    sf_solver_reg(dsrtomo_oper,sf_cgstep,sf_igrad2_lop,2*nw,nw,nt,dw,dt,cgiter,eps,"verb",verb,"end");
		    sf_cgstep_close();
		}

		if (grad != NULL) {
		    if (velocity) {
			for (iw=0; iw < nw; iw++)
			    dv[iw] = -dw[iw]/(2.*sqrtf(w[iw])*(w[iw]+dw[iw]/2.));
			
			sf_floatwrite(dv,nw,grad);
		    } else {
			sf_floatwrite(dw,nw,grad);
		    }
		}

		/* line search */
		gama = 0.5;
		for (count=0; count < 5; count++) {
		    
		    /* update slowness */
		    for (iw=0; iw < nw; iw++)
			w1[iw] = (w[iw]+gama*dw[iw])*(w[iw]+gama*dw[iw])/w[iw];
		    
		    /* compute new misfit */
		    dsreiko_fastmarch(t,w1,f,th);
		    dsreiko_mirror(t);
		    
		    for (it=0; it < nt; it++) {
			if (m0 == NULL || m0[it] == 1)
			    dt[it] = t0[it]-t[it];
			else
			    dt[it] = 0.;
		    }
		    
		    rhsnorm = cblas_snrm2(nt,dt,1);
		    rate = rhsnorm/rhsnorm1;
		    
		    if (rate < 1.) {
			for (iw=0; iw < nw; iw++)
			    w[iw] = w1[iw];
			
			rhsnorm1 = rhsnorm;
			rate = rhsnorm1/rhsnorm0;
			break;
		    }
		    
		    gama *= 0.5;
		}
		
		if (count == 5) {
		    sf_warning("Line-search failure at iteration %d of %d.",iter+1,niter);
		    break;
		}
		
		sf_warning("L2 misfit after iteration %d of %d: %g (line-search %d)",iter+1,niter,rate,count);
	    }
	    
	    /* convert to velocity */
	    if (velocity) {
		for (iw=0; iw < nw; iw++) {
		    w[iw] = 1./sqrtf(w[iw]);
		}
	    }
	    
	    sf_floatwrite(w,nw,out);
	    
	    break;
    }
    
    exit(0);
}
