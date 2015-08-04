/* First-arrival Traveltime Tomography */
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

#include "fastmarch.h"
#include "fatomo.h"

#include "l1.h"
#include "l1step.c"

int main(int argc, char* argv[])
{
    bool adj, velocity, l1norm, plane[3], verb;
    int dim, i, count, n[SF_MAX_DIM], it, nt, **m, nrhs, is, nshot=1, *flag, order, iter, niter, stiter, *k, nfreq, nmem;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], **t, *t0, *s, *temps, **source, *rhs, *ds;
    float rhsnorm, rhsnorm0, rhsnorm1, rate, eps, gama;
    char key[4], *what;
    sf_file sinp, sout, shot, time, reco, rece, topo, grad, norm;
    sf_weight weight=NULL;

    sf_init(argc,argv);
    sinp = sf_input("in");
    sout = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="tomo";
    /* what to compute (default tomography) */

    switch (what[0]) {
	case 'l': /* linear operator */

	    if (NULL == sf_getstring("time"))
		sf_error("Need time=");
	    time = sf_input("time");

	    /* read operator dimension from time table */
	    dim = sf_filedims(time,n);
	    
	    nt = 1;
	    for (i=0; i < 3; i++) {
		sprintf(key,"d%d",i+1);
		if (!sf_histfloat(time,key,d+i)) sf_error("No %s= in input",key);
		sprintf(key,"o%d",i+1);
		if (!sf_histfloat(time,key,o+i)) o[i]=0.;
		nt *= n[i]; plane[i] = false;
	    }
	    if (dim < 3) {
		n[2] = 1; o[2] = o[1]; d[2] = d[1]; plane[2] = false;
	    }

	    dim = 2;

	    /* read in shot file */
	    if (NULL == sf_getstring("shot"))
		sf_error("Need source shot=");
	    shot = sf_input("shot");
	    
	    if (!sf_histint(shot,"n2",&nshot)) nshot=1;
	    sf_fileclose(shot);

	    /* read in receiver file */
	    m = sf_intalloc2(nt,nshot);
	    if (NULL == sf_getstring("receiver")) {
		for (is=0; is < nshot; is++) {
		    for (it=0; it < nt; it++) {
			m[is][it] = 1;
		    }
		}
	    } else {
		rece = sf_input("receiver");
		sf_intread(m[0],nt*nshot,rece);
		sf_fileclose(rece);
	    }

	    /* number of right-hand side */
	    nrhs = 0;
	    for (is=0; is < nshot; is++) {
		for (it=0; it < nt; it++) {
		    if (m[is][it] == 1) nrhs++;
		}
	    }
	    rhs = sf_floatalloc(nrhs);

	    t = sf_floatalloc2(nt,nshot);
	    sf_floatread(t[0],nt*nshot,time);

	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */

	    /* initialize fatomo */
	    fatomo_init(dim,n,d,nshot);

	    /* set operators */
	    fatomo_set(t,m);

	    t0 = sf_floatalloc(nt);

	    if (adj) {
		sf_floatread(rhs,nrhs,sinp);

		fatomo_lop(true,false,nt,nrhs,t0,rhs);

		sf_putint(sout,"n1",nt);
		sf_putint(sout,"n2",1);
		sf_putint(sout,"n3",1);
		sf_floatwrite(t0,nt,sout);
	    } else {
		sf_floatread(t0,nt,sinp);

		fatomo_lop(false,false,nt,nrhs,t0,rhs);
		
		sf_putint(sout,"n1",nrhs);
		sf_putint(sout,"n2",1);
		sf_putint(sout,"n3",1);
		sf_floatwrite(rhs,nrhs,sout);
	    }

	    break;
	    
	case 't': /* tomography */

	    /* read input dimension */
	    dim = sf_filedims(sinp,n);
	    
	    nt = 1;
	    for (i=0; i < dim; i++) {
		sprintf(key,"d%d",i+1);
		if (!sf_histfloat(sinp,key,d+i)) sf_error("No %s= in input",key);
		sprintf(key,"o%d",i+1);
		if (!sf_histfloat(sinp,key,o+i)) o[i]=0.;
		nt *= n[i]; plane[i] = false;
	    }
	    if (dim < 3) {
		n[2] = 1; o[2] = o[1]; d[2] = d[1]; plane[2] = false;
	    }
	    
	    /* read initial guess */
	    s = sf_floatalloc(nt);
	    sf_floatread(s,nt,sinp);
	    
	    if (!sf_getbool("velocity",&velocity)) velocity=true;
	    /* if y, the input is velocity; n, slowness squared */

	    if (velocity) {
		for (it=0; it < nt; it++) {
		    s[it] = 1./s[it]*1./s[it];
		}
	    }

	    /* allocate memory for temporary data */
	    ds    = sf_floatalloc(nt);
	    flag  = sf_intalloc(nt);
	    
	    temps = sf_floatalloc(nt);
	    for (it=0; it < nt; it++) {
		temps[it] = s[it];
	    }
	    
	    if (!sf_getbool("l1norm",&l1norm)) l1norm=false;
	    /* norm for minimization (default L2 norm) */

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

	    /* allocate memory for time table */
	    t = sf_floatalloc2(nt,nshot);

	    /* read in receiver file */
	    m = sf_intalloc2(nt,nshot);
	    if (NULL == sf_getstring("receiver")) {
		for (is=0; is < nshot; is++) {
		    for (it=0; it < nt; it++) {
			m[is][it] = 1;
		    }
		}
	    } else {
		rece = sf_input("receiver");
		sf_intread(m[0],nt*nshot,rece);
		sf_fileclose(rece);
	    }
	    
	    /* number of right-hand side */
	    nrhs = 0;
	    for (is=0; is < nshot; is++) {
		for (it=0; it < nt; it++) {
		    if (m[is][it] == 1) nrhs++;
		}
	    }
	    rhs = sf_floatalloc(nrhs);
	    
	    /* read in record file */
	    if (NULL == sf_getstring("record"))
		sf_error("Need data record=");
	    reco = sf_input("record");

	    t0 = sf_floatalloc(nrhs);
	    sf_floatread(t0,nrhs,reco);
	    sf_fileclose(reco);

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

	    if (!sf_getint("niter",&niter)) niter=10;
	    /* number of slowness inversion iterations */

	    if (!sf_getint("stiter",&stiter)) stiter=200;
	    /* number of step iterations */

	    if (!sf_getfloat("eps",&eps)) eps=0.;
	    /* regularization parameter */

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

	    /* initialize fatomo */
	    fatomo_init(dim,n,d,nshot);

	    /* initialize 2D gradient operator */
	    sf_igrad2_init(n[0],n[1]);

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

	    /* initial misfit */
	    fastmarch_init(n[2],n[1],n[0]);
	    
	    i = 0;
	    for (is=0; is < nshot; is++) {
		fastmarch(t[is],s,flag,plane,
			  n[2],n[1],n[0],o[2],o[1],o[0],d[2],d[1],d[0],
			  source[is][2],source[is][1],source[is][0],1,1,1,order);

		for (it=0; it < nt; it++) {
		    if (m[is][it] == 1) {
			rhs[i] = t0[i]-t[is][it];
			i++;
		    }
		}
	    }
	    
	    fastmarch_close();
	    
	    /* calculate L2 data-misfit */
	    rhsnorm0 = cblas_snrm2(nrhs,rhs,1);
	    rhsnorm = rhsnorm0;
	    rhsnorm1 = rhsnorm;
	    rate = rhsnorm1/rhsnorm0;

	    if (l1norm)
		sf_warning("L1 misfit after iteration 0 of %d: %g",niter,rate);
	    else
		sf_warning("L2 misfit after iteration 0 of %d: %g",niter,rate);

	    if (norm != NULL) sf_floatwrite(&rate,1,norm);

	    /* iterations over inversion */
	    for (iter=0; iter < niter; iter++) {
		
		/* clean-up */
		for (it=0; it < nt; it++) {
		    ds[it] = 0.;
		}

		/* prepare for CG */
		fatomo_set(t,m);

		/* solve ds */
		if (l1norm) {
		    /*
		      sf_solver_reg(fatomo_lop,l1step,sf_igrad2_lop,2*nt, nt,nrhs,ds,rhs,stiter,eps,"verb",verb,"end");
		    */
		    sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt,nt,nrhs,ds,rhs,stiter,eps,"wght",weight,"nfreq",nfreq,"nmem",nmem,"verb",verb,"end");
		    
		    /*
		      l1step_close();
		    */
		    
		    sf_cgstep_close();
		} else {
		    sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt,nt,nrhs,ds,rhs,stiter,eps,"verb",verb,"end");
			
		    sf_cgstep_close();
		}
		
		/* line search */
		gama = 1.;
		for (count=0; count < 10; count++) {
		    
		    /* update slowness */
		    for (it=0; it < nt; it++) {
			if (k == NULL || k[it] != 1)
			    temps[it] = (s[it]+gama*ds[it])*(s[it]+gama*ds[it])/s[it];
		    }

		    /* forward fast-marching for stencil time */
		    fastmarch_init(n[2],n[1],n[0]);
		    
		    i = 0;
		    for (is=0; is < nshot; is++) {
			fastmarch(t[is],temps,flag,plane,
				  n[2],n[1],n[0],o[2],o[1],o[0],d[2],d[1],d[0],
				  source[is][2],source[is][1],source[is][0],1,1,1,order);
			
			for (it=0; it < nt; it++) {
			    if (m[is][it] == 1) {
				rhs[i] = t0[i]-t[is][it];
				i++;
			    }
			}
		    }
		    
		    fastmarch_close();
		    
		    rhsnorm = cblas_snrm2(nrhs,rhs,1);
		    rate = rhsnorm/rhsnorm1;
		    
		    if (rate < 1.) {
			for (it=0; it < nt; it++) {
			    s[it] = temps[it];
			}
			rhsnorm1 = rhsnorm;
			rate = rhsnorm1/rhsnorm0;
			break;
		    }
		    
		    gama *= 0.5;
		}
		
		if (count == 10) {
		    sf_warning("Line-search Failure. Iteration terminated at %d of %d.",iter+1,niter);
		    sf_warning("Dimensions for GRAD and NORM need to be fixed before read.");
		    break;
		}

		if (l1norm)
		    sf_warning("L1 misfit after iteration %d of %d: %g (line-search %d)",iter+1,niter,rate,count);
		else
		    sf_warning("L2 misfit after iteration %d of %d: %g (line-search %d)",iter+1,niter,rate,count);
		
		if (grad != NULL) sf_floatwrite(ds,nt,grad);
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
