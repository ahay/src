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
#include "trianglen.h"

int main(int argc, char* argv[])
{
    bool velocity, shape, verb, weight;
    int dim, i, n[SF_MAX_DIM], rect[SF_MAX_DIM], it, nt, **m, is, nshot, order, seg;
    int iter, niter, stiter, istep, nstep, *k, nrhs, **rhslist, nrecv, *mp;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], **t, **t0, *s, *temps, *dv=NULL, **source, *rhs, *ds, *p=NULL, **modl, **ray, *wght=NULL;
    float tol, rhsnorm, rhsnorm0, rhsnorm1, rate, eps, step, pow;
    char key[6], *what;
    sf_file sinp, sout, shot, reco, recv, topo, grad, norm, rayd, time, prec;
    
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

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */
    
    /* convert to slowness squared */
    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./s[it]*1./s[it];
	}
	
	dv = sf_floatalloc(nt);
    }
    
    /* allocate memory for temporary data */
    ds = sf_floatalloc(nt);
    
    temps = sf_floatalloc(nt);
    for (it=0; it < nt; it++) {
	temps[it] = s[it];
    }

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
    sf_fileclose(shot);
  
    /* read in receiver list */
    if (NULL == sf_getstring("recv"))
	sf_error("Need receiver list recv=");
    recv = sf_input("recv");

    if (!sf_histint(recv,"n1",&nrecv)) 
	sf_error("No nrecv in receiver list.");

    m = sf_intalloc2(nrecv,nshot);
    sf_intread(m[0],nrecv*nshot,recv);
    sf_fileclose(recv);

    /* number of right-hand side */
    rhslist = sf_intalloc2(2,nshot);
    
    nrhs = 0;
    for (is=0; is < nshot; is++) {
	/* rhslist[is][0]: where this shot starts in rhs vector */
	rhslist[is][0] = nrhs;
	
	for (it=0; it < nrecv; it++) {
	    if (m[is][it] >= 0)
		nrhs++;
	}
	
	/* rhslist[is][1]: how many receivers there are for this shot */
	rhslist[is][1] = (is==0)? nrhs: (nrhs-rhslist[is-1][0]-rhslist[is-1][1]);
    }
    rhs = sf_floatalloc(nrhs);

    /* read in record list */
    if (NULL == sf_getstring("reco"))
	sf_error("Need record list reco=");
    reco = sf_input("reco");
    
    t0 = sf_floatalloc2(nrecv,nshot);
    sf_floatread(t0[0],nrecv*nshot,reco);
    sf_fileclose(reco);

    /* read in topography file */
    if (NULL != sf_getstring("topo")) {
	topo = sf_input("topo");
	k = sf_intalloc(nt);
	sf_intread(k,nt,topo);
	sf_fileclose(topo);
    } else {
	topo = NULL;
	k = NULL;
    }
    
    /* read model mask */
    if (NULL == sf_getstring("prec")) {
	prec = NULL;
	mp = NULL;
    } else {
	prec = sf_input("prec");
	mp = sf_intalloc(nt);
	sf_intread(mp,nt,prec);
	sf_fileclose(prec);
    }

    if (!sf_getint("order",&order)) order=2;
    /* fast marching accuracy order */
    
    if (!sf_getint("seg",&seg)) seg=3;
    /* maximum number of segments of topography */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of slowness inversion iterations */
    
    if (!sf_getint("stiter",&stiter)) stiter=100;
    /* number of inner CG iterations (for both Ticknov and Shaping) */
    
    if (!sf_getint("nstep",&nstep)) nstep=10;
    /* number of linesearch */

    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameter (for both Ticknov and Shaping) */
    
    if (!sf_getbool("weight",&weight)) weight=false;
    /* data weighting */

    if (weight) {
	wght = sf_floatalloc(nrhs);

	if (!sf_getfloat("pow",&pow)) pow=2.;
	/* power raised for data weighting */
    }

    /* output gradient at each iteration */
    if (NULL != sf_getstring("grad")) {
	grad = sf_output("grad");
	sf_putint(grad,"n3",n[2]);
	sf_putfloat(grad,"d3",d[2]);
	sf_putfloat(grad,"o3",o[2]);
	sf_putint(grad,"n4",niter);
    } else {
	grad = NULL;
    }
    
    /* output ray density/coverage at each iteration */
    if (NULL != sf_getstring("rayd")) {
	rayd = sf_output("rayd");
	sf_putint(rayd,"n3",n[2]);
	sf_putfloat(rayd,"d3",d[2]);
	sf_putfloat(rayd,"o3",o[2]);
	sf_putint(rayd,"n4",nshot);
	sf_putint(rayd,"n5",niter+1);
	ray = sf_floatalloc2(nt,nshot);
    } else {
	rayd = NULL;
	ray = NULL;
    }

    /* output forward-modeled record at each iteration */
    if (NULL != sf_getstring("time")) {
	time = sf_output("time");
	sf_putint(time,"n1",nrecv);
	sf_putint(time,"n2",nshot);
	sf_putint(time,"n3",niter+1);
	modl = sf_floatalloc2(nrecv,nshot);
    } else {
	time = NULL;
	modl = NULL;
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
	if (topo != NULL) {
	    trianglen_init(dim,rect,n);
	    trianglen_topo(k,seg);
	    sf_repeat_init(nt,1,trianglen_lop);
	} else {
	    sf_trianglen_init(dim,rect,n);
	    sf_repeat_init(nt,1,sf_trianglen_lop);
	}

	sf_conjgrad_init(nt,nt,nrhs,nrhs,eps,tol,verb,false);
	p = sf_floatalloc(nt);
    }
    
    /* initialize fatomo */
    fatomo_init(dim,n,o,d,order,nshot,rhslist,m,t0,wght,mp);

    /* initial misfit */
    fatomo_fastmarch(s,t,source,rhs);

    /* data weighting */
    if (weight) {
	for (i=0; i < nrhs; i++) {
	    wght[i] = 1.;
	}
    }
    
    /* output forward-modeled record */
    if (time != NULL) {
	for (is=0; is < nshot; is++) {
	    for (it=0; it < nrecv; it++) {
		if (it < rhslist[is][1])
		    modl[is][it] = t0[is][it]-rhs[rhslist[is][0]+it];
		else
		    modl[is][it] = 0.;
	    }
	}

	sf_floatwrite(modl[0],nrecv*nshot,time);
    }
    
    /* output ray density/coverage */
    if (rayd != NULL) {
	fatomo_ray(ray);
	sf_floatwrite(ray[0],nt*nshot,rayd);
    }

    /* calculate L2 data-misfit */
    rhsnorm0 = cblas_snrm2(nrhs,rhs,1);
    rhsnorm = rhsnorm0;
    rhsnorm1 = rhsnorm;
    rate = rhsnorm1/rhsnorm0;
    
    sf_warning("L2 misfit after iteration 0 of %d:\t %g",niter,rate);
    
    if (norm != NULL) sf_floatwrite(&rate,1,norm);
    
    switch (what[0]) {
	case 'l': /* linear operator */

	    fatomo_lop(true,false,nt,nrhs,ds,rhs);

	    if (velocity) {
		for (it=0; it < nt; it++) {		    
		    dv[it] = -ds[it]/(sqrtf(s[it])*(s[it]+ds[it]));
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
		
		/* clean-up */
		for (it=0; it < nt; it++) {
		    ds[it] = 0.;
		}

		/* solve ds */
		if (shape) {
		    sf_conjgrad(NULL,fatomo_lop,sf_repeat_lop,p,ds,rhs,stiter);
		} else {
		    sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt,nt,nrhs,ds,rhs,stiter,eps,"verb",verb,"end");
		    sf_cgstep_close();
		}
		
		/* output computed gradient (before line-search) */
		if (grad != NULL) {
		    if (velocity) {
			for (it=0; it < nt; it++) {			    
			    dv[it] = -ds[it]/(sqrtf(s[it])*(s[it]+ds[it]));
			}
			sf_floatwrite(dv,nt,grad);
		    } else {
			sf_floatwrite(ds,nt,grad);
		    }
		}

		/* line search */
		for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {
		    
		    /* update slowness */
		    for (it=0; it < nt; it++) {
			if (k == NULL || k[it] == 1)
			    temps[it] = (s[it]+step*ds[it])*(s[it]+step*ds[it])/s[it];
		    }
		    
		    /* forward fast-marching for stencil time */		    
		    fatomo_fastmarch(temps,t,source,rhs);

		    if (weight) {
			for (i=0; i < nrhs; i++) {
			    wght[i] = 1./expf(fabsf(powf(iter+1.,pow)*rhs[i]));
			    rhs[i] *= wght[i];
			}
		    }

		    rhsnorm = cblas_snrm2(nrhs,rhs,1);
		    rate = rhsnorm/rhsnorm1;
		    
		    /* break */
		    if (rate < 1.) {
			for (it=0; it < nt; it++) {
			    s[it] = temps[it];
			}
			rhsnorm1 = rhsnorm;
			rate = rhsnorm1/rhsnorm0;

			break;
		    }
		}
		
		/* output forward-modeled record */
		if (time != NULL) {
		    for (is=0; is < nshot; is++) {
			for (it=0; it < nrecv; it++) {
			    if (it < rhslist[is][1])
				modl[is][it] = t0[is][it]-rhs[rhslist[is][0]+it];
			    else
				modl[is][it] = 0.;
			}
		    }
		    
		    sf_floatwrite(modl[0],nrecv*nshot,time);
		}
		
		/* output ray density/coverage */
		if (rayd != NULL) {
		    fatomo_ray(ray);
		    sf_floatwrite(ray[0],nt*nshot,rayd);
		}

		if (istep == 10) {
		    sf_warning("Line-search Failure. Iteration terminated at %d of %d.",iter+1,niter);
		    sf_warning("Dimensions for NORM need to be fixed before read.");
		    break;
		}

		sf_warning("L2 misfit after iteration %d of %d:\t %g (istep = %d).",iter+1,niter,rate,istep);

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
