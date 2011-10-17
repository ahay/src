/* First-arrival Traveltime Tomography (OMP) */
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

#include "upgradomp.h"
#include "fatomoomp.h"

#include "l1.h"
#include "l1step.c"

int main(int argc, char* argv[])
{
    bool velocity, l1norm, shape, verb;
    int dim, i, n[SF_MAX_DIM], rect[SF_MAX_DIM], it, nt, **m, is, nshot, order;
    int iter, niter, stiter, istep, nstep, *k, nfreq, nmem, nrhs, **rhslist, nrecv;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], **t, **t0, *s, *temps, *dv=NULL, **source, *rhs, *ds, *p=NULL;
    float tol, rhsnorm, rhsnorm0, rhsnorm1, rate, eps, step, min, max, gama, ratio;
    char key[6], *what;
    sf_file sinp, sout, shot, reco, recv, topo, grad, norm, steep, titer;
    sf_weight weight=NULL;

    sf_init(argc,argv);
    sinp = sf_input("in");
    sout = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="tomo";
    /* what to compute (default tomography) */

    /* read input dimension */
    dim = sf_filedims(sinp,n);
    
    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(sinp,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(sinp,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
    }
    
    /* read initial guess */
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,sinp);

    if (!sf_getfloat("min",&min)) min=0.;
    /* minimum allowed velocity */

    if (!sf_getfloat("max",&max)) max=SF_HUGE;
    /* maximum allowed velocity */

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */
    
    /* convert to slowness squared */
    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./s[it]*1./s[it];
	}
	
	dv = sf_floatalloc(nt);
	
	min = (min==0.)?SF_HUGE:(1./min);
	max = (max==SF_HUGE)?0.:(1./max);
    }
    
    /* allocate memory for temporary data */
    ds = sf_floatalloc(nt);
    
    temps = sf_floatalloc(nt);
    for (it=0; it < nt; it++) {
	temps[it] = s[it];
    }

    if (!sf_getbool("l1norm",&l1norm)) l1norm=false;
    /* norm for minimization (default L2 norm) */
    
    if (!sf_getbool("shape",&shape)) shape=false;
    /* regularization (default Tikhnov) */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    /* read in shot file */
    if (NULL == sf_getstring("shot"))
	sf_error("Need source shot=");
    shot = sf_input("shot");
    
    if (!sf_histint(shot,"n2",&nshot)) nshot=1;
    
    source = sf_floatalloc2(3,nshot);
    sf_floatread(source[0],3*nshot,shot);
/*
    sf_fileclose(shot);
*/  
  
    /* read in receiver list */
    if (NULL == sf_getstring("receiver"))
	sf_error("Need list receiver=");
    recv = sf_input("receiver");

    if (!sf_histint(recv,"n1",&nrecv)) 
	sf_error("No nrecv in receiver list.");

    m = sf_intalloc2(nrecv,nshot);
    sf_intread(m[0],nrecv*nshot,recv);
/*
    sf_fileclose(recv);
*/

    /* number of right-hand side */
    rhslist = sf_intalloc2(2,nshot);
    
    nrhs = 0;
    for (is=0; is < nshot; is++) {
	rhslist[is][1] = nrhs;
	
	for (it=0; it < nrecv; it++) {
	    if (m[is][it] >= 0)
		nrhs++;
	}
	
	rhslist[is][2] = (is==0)?nrhs
	    :rhslist[is][1]-rhslist[is-1][1];
    }
    rhs = sf_floatalloc(nrhs);

    /* read in record list */
    if (NULL == sf_getstring("record"))
	sf_error("Need list record=");
    reco = sf_input("record");
    
    t0 = sf_floatalloc2(nrecv,nshot);
    sf_floatread(t0[0],nrecv*nshot,reco);
/*
    sf_fileclose(reco);
*/    

    /* read in topography file */
    if (NULL != sf_getstring("topo")) {
	topo = sf_input("topo");
	k = sf_intalloc(nt);
	sf_intread(k,nt,topo);
	sf_fileclose(topo);
    } else {
	k = NULL;
    }
    
    if (!sf_getint("order",&order)) order=2;
    /* fast marching accuracy order */
    
    if (!sf_getint("niter",&niter)) niter=1;
    /* number of slowness inversion iterations */
    
    if (!sf_getint("stiter",&stiter)) stiter=100;
    /* number of inner CG iterations (for both Ticknov and Shaping) */
    
    if (!sf_getint("nstep",&nstep)) nstep=10;
    /* number of linesearch */

    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameter (for both Ticknov and Shaping) */
    
    /* output gradient at each iteration */
    if (NULL != sf_getstring("gradient")) {
	grad = sf_output("gradient");
	sf_putint(grad,"n3",n[2]);
	sf_putfloat(grad,"d3",d[2]);
	sf_putfloat(grad,"o3",o[2]);
	sf_putint(grad,"n4",niter);
    } else {
	grad = NULL;
    }
    
    /* output steepest-descent at each iteration */
    if (NULL != sf_getstring("steep")) {
	steep = sf_output("steep");
	sf_putint(steep,"n3",n[2]);
	sf_putfloat(steep,"d3",d[2]);
	sf_putfloat(steep,"o3",o[2]);
	sf_putint(steep,"n4",niter);
    } else {
	steep = NULL;
    }

    /* output time at each iteration */
    if (NULL != sf_getstring("titer")) {
	titer = sf_output("titer");
	sf_putint(titer,"n3",n[2]);
	sf_putfloat(titer,"d3",d[2]);
	sf_putfloat(titer,"o3",o[2]);
	sf_putint(titer,"n4",nshot);
	sf_putint(titer,"n5",niter);
    } else {
	titer = NULL;
    }

    /* output misfit L2 norm at each iteration */
    if (NULL != sf_getstring("misnorm")) {
	norm = sf_output("misnorm");
	sf_putint(norm,"n1",niter+1);
	sf_putfloat(norm,"d1",1.);
	sf_putfloat(norm,"o1",0.);
	sf_putint(norm,"n2",1);
	sf_putint(norm,"n3",1);
    } else {
	norm = NULL;
    }
    
    /* allocate memory for time table */
    t = sf_floatalloc2(nt,nshot);
    
    /* initialize 2D gradient operator */
    sf_igrad2_init(n[0],n[1]);
    
    /* initialize L1-norm operator */
    if (l1norm) {
	/*
	  if (!sf_getfloat("perc",&perc)) perc=90.;
	  
	  l1_init(nt,stiter,perc,false);
	*/
	if (!sf_getint("nfreq",&nfreq)) nfreq=1;
	/* l1-norm weighting nfreq */
	
	if (!sf_getint("nmem",&nmem)) nmem=1;
	/* l1-norm weighting nmem */
	
	weight = sf_l1;
	sf_irls_init(nt);
    }
    
    /* initialize shaping operator */
    if (shape) {
	if (!sf_getfloat("tol",&tol)) tol=1.e-6;
	/* tolerance for shaping regularization */

	for (i=0; i < dim; i++) {
	    sprintf(key,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	    /*( rect#=(1,1,...) smoothing radius on #-th axis )*/
	}

	/* triangle smoothing operator */
	sf_trianglen_init(dim,rect,n);
	sf_repeat_init(nt,1,sf_trianglen_lop);

	sf_conjgrad_init(nt,nt,nrhs,nrhs,eps,tol,verb,false);
	p = sf_floatalloc(nt);
    }
    
    /* initialize fatomo */
    fatomo_init(dim,n,o,d,order,nshot,rhslist,m,t0);

    /* initial misfit */
    fatomo_fastmarch(s,t,source,rhs);
    
    if (titer != NULL)
	sf_floatwrite(t[0],nt*nshot,titer);

    /* calculate L2 data-misfit */
    rhsnorm0 = cblas_snrm2(nrhs,rhs,1);
    rhsnorm = rhsnorm0;
    rhsnorm1 = rhsnorm;
    rate = rhsnorm1/rhsnorm0;
    
    if (l1norm)
	sf_warning("L1 misfit after iteration 0 of %d:\t %g",niter,rate);
    else
	sf_warning("L2 misfit after iteration 0 of %d:\t %g",niter,rate);
    
    if (norm != NULL) sf_floatwrite(&rate,1,norm);
    
    switch (what[0]) {
	case 'l': /* linear operator */

	    fatomo_lop(true,false,nt,nrhs,ds,rhs);

	    if (velocity) {
		for (it=0; it < nt; it++) {
		    dv[it] = -ds[it]/(s[it]+sqrtf(s[it])*ds[it]);
		}
		sf_floatwrite(dv,nt,sout);
	    } else {
		sf_floatwrite(ds,nt,sout);
	    }
	    
	    break;
	    
	case 't': /* tomography */
	    
	    /* iterations over inversion */
	    for (iter=0; iter < niter; iter++) {
		
		/* clean-up */
		for (it=0; it < nt; it++) {
		    ds[it] = 0.;
		}
		
		if (steep != NULL) {
		    fatomo_lop(true,false,nt,nrhs,ds,rhs);

		    if (velocity) {
			for (it=0; it < nt; it++) {
			    dv[it] = -ds[it]/(s[it]+sqrtf(s[it])*ds[it]);
			}
			sf_floatwrite(dv,nt,steep);
		    } else {
			sf_floatwrite(ds,nt,steep);
		    }
		}

		/* clean-up */
		for (it=0; it < nt; it++) {
		    ds[it] = 0.;
		}

		/* solve ds */
		if (l1norm) {
		    /* sf_solver_reg(fatomo_lop,l1step,sf_igrad2_lop,2*nt, nt,nrhs,ds,rhs,stiter,eps,"verb",verb,"end"); */
		    sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt,nt,nrhs,ds,rhs,stiter,eps,"wght",weight,"nfreq",nfreq,"nmem",nmem,"verb",verb,"end");
		    
		    /* l1step_close(); */
		    sf_cgstep_close();
		} else {
		    if (shape) {
			sf_conjgrad(NULL,fatomo_lop,sf_repeat_lop,p,ds,rhs,stiter);
		    } else {
			sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt,nt,nrhs,ds,rhs,stiter,eps,"verb",verb,"end");
			sf_cgstep_close();
		    }
		}

		/* convert to slowness perturbation */
		for (it=0; it < nt; it++) {
		    ds[it] = ds[it]/sqrtf(s[it]);
		}
		
		if (grad != NULL) {
		    if (velocity) {
			for (it=0; it < nt; it++) {
			    dv[it] = -ds[it]/(s[it]+sqrtf(s[it])*ds[it]);
			}
			sf_floatwrite(dv,nt,grad);
		    } else {
			sf_floatwrite(ds,nt,grad);
		    }
		}

		/* avoid overshot */
		gama = 1.;

		for (it=0; it < nt; it++) {
		    /* go below minimum allowed value */
		    if (ds[it]>0. && (sqrtf(s[it])+ds[it])>min) {
			ratio = (min-sqrtf(s[it]))/ds[it];
			gama = (gama<ratio)?gama:ratio;
		    }

		    /* go beyond maximum allowed value */
		    if (ds[it]<0. && (sqrtf(s[it])+ds[it])<max) {
			ratio = (max-sqrtf(s[it]))/ds[it];
			gama = (gama<ratio)?gama:ratio;
		    }
		}

		/* line search */
		for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {
		    
		    /* update slowness */
		    for (it=0; it < nt; it++) {
			if (k == NULL || k[it] != 1)
			    temps[it] = (sqrtf(s[it])+gama*step*ds[it])*(sqrtf(s[it])+gama*step*ds[it]);
		    }
		    
		    /* forward fast-marching for stencil time */		    
		    fatomo_fastmarch(temps,t,source,rhs);

		    rhsnorm = cblas_snrm2(nrhs,rhs,1);
		    rate = rhsnorm/rhsnorm1;
		    
		    /* break */
		    if (rate < 1.) {
			for (it=0; it < nt; it++) {
			    s[it] = temps[it];
			}
			rhsnorm1 = rhsnorm;
			rate = rhsnorm1/rhsnorm0;

			if (titer != NULL)
			    sf_floatwrite(t[0],nt*nshot,titer);

			break;
		    }
		}
		
		if (istep == 10) {
		    sf_warning("Line-search Failure. Iteration terminated at %d of %d.",iter+1,niter);
		    sf_warning("Dimensions for NORM need to be fixed before read.");
		    break;
		}

		if (l1norm)
		    sf_warning("L1 misfit after iteration %d of %d:\t %g (gama = %g, istep = %d).",iter+1,niter,rate,gama,istep);
		else
		    sf_warning("L2 misfit after iteration %d of %d:\t %g (gama = %g, istep = %d).",iter+1,niter,rate,gama,istep);

		if (norm != NULL) sf_floatwrite(&rate,1,norm);
	    }
	    
	    /* convert to velocity */
	    if (velocity) {
		for (it=0; it < nt; it++) {
		    s[it] = 1./sqrtf(s[it]);
		}
	    }
	    
	    sf_floatwrite(s,nt,sout);
	    
       	    break;
    }
    
    exit(0);
}
