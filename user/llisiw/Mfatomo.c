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

#include "upgrad.h"
#include "fastmarch.h"
#include "fatomo.h"

int main(int argc, char* argv[])
{
    bool adj, velocity, plane[3];
    int dim, i, n[SF_MAX_DIM], it, nt, **m, nrhs, is, nshot, *flag, order, iter, niter, cgiter, *k;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *t, **t0, *s, **source, *rhs, *ds, *gs, air, *x0;
    float rhsnorm0, rhsnorm, rate, eps;
    char key[4], *what;
    upgrad upg;
    sf_file sinp, sout, shot, time, reco, rece, topo, grad, norm;

    sf_init(argc,argv);
    sinp = sf_input("in");
    sout = sf_output("out");

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

    if (NULL == (what = sf_getstring("what"))) what="tomo";
    /* what to compute (default tomography) */

    s = sf_floatalloc(nt);
    t = sf_floatalloc(nt);
 
    sf_floatread(s,nt,sinp);

    switch (what[0]) {
	case 'l': /* linear operator */
	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */

	    if (NULL == sf_getstring("time"))
		sf_error("Need background time=");
	    time = sf_input("time");
	    sf_floatread(t,nt,time);
	    sf_fileclose(time);

	    /* set stencil */
	    upg = upgrad_init(dim,n,d);
	    upgrad_set(upg,t);

	    if (adj) {
		upgrad_inverse(upg,t,s,NULL);
	    } else {
		upgrad_solve(upg,s,t,NULL);
	    }

	    sf_floatwrite(t,nt,sout);

	    break;
	    
	case 't': /* tomography */
	    ds   = sf_floatalloc(nt);
	    gs   = sf_floatalloc(nt);
	    flag = sf_intalloc(nt);
	    
	    /* read in shot file */
	    if (NULL == sf_getstring("shot"))
		sf_error("Need source shot=");
	    shot = sf_input("shot");

	    if (!sf_histint(shot,"n2",&nshot)) nshot=1;

	    source = sf_floatalloc2(3,nshot);
	    sf_floatread(source[0],3*nshot,shot);
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

	    nrhs = 0;
	    for (it=0; it < nt; it++) {
		if (m[0][it] == 1) nrhs++;
	    }
	    rhs = sf_floatalloc(nrhs);

	    /* read in record file */
	    if (NULL == sf_getstring("record"))
		sf_error("Need data record=");
	    reco = sf_input("record");

	    /* read in topography file */
	    if (NULL != sf_getstring("topo")) {
		topo = sf_input("topo");
		k = sf_intalloc(nt);
		sf_intread(k,nt,topo);
		sf_fileclose(topo);

		x0 = sf_floatalloc(nt);

		if (!sf_getfloat("air",&air)) air = 0.01;
		/* air velocity (for fixed topo inversion) */

		for (it=0; it < nt; it++) {
		    if (k[it] == 0) {
			k[it] = 1;
			x0[it] = air;
		    } else {
			k[it] = 0;
			x0[it] = 0.;
		    }
		}
		    
	    }

	    t0 = sf_floatalloc2(nrhs,nshot);
	    sf_floatread(t0[0],nrhs*nshot,reco);
	    sf_fileclose(reco);

	    if (!sf_getint("order",&order)) order=2;
	    /* fast marching accuracy order */

	    if (!sf_getbool("velocity",&velocity)) velocity=true;
	    /* if y, the input is velocity; n, slowness squared */

	    if (!sf_getint("niter",&niter)) niter=10;
	    /* number of slowness inversion iterations */

	    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
	    /* number of conjugate gradient iterations */

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

	    if (velocity) {
		for (it=0; it < nt; it++) {
		    s[it] = 1./s[it]*1./s[it];
		}
	    }

	    rhsnorm0 = 0.;
	    fatomo_init(dim,n,d);

	    if (!sf_getfloat("eps",&eps)) eps=0.;
	    /* regularization parameter */

	    sf_igrad2_init(n[0],n[1]);

	    /* iterations over inversion */
	    for (iter=0; iter < niter; iter++) {

		/* clean-up */
		for (it=0; it < nt; it++) {
		    gs[it] = 0.;
		}

		rhsnorm = 0.;

		/* loop over all shots */
		for (is=0; is < nshot; is++) {
		    sf_warning("shot %d of %d, iteration %d of %d",is+1,nshot,iter+1,niter);

		    for (it=0; it < nt; it++) {
			t[it] = 0.;
			ds[it] = 0.;
		    }
		    
		    /* forward fast-marching for stencil time */
		    fastmarch_init(n[2],n[1],n[0]);

		    fastmarch(t,s,flag,plane,
			      n[2],n[1],n[0],o[2],o[1],o[0],d[2],d[1],d[0],
			      source[is][2],source[is][1],source[is][0],1,1,1,order);

		    fastmarch_close();
/*
		    sf_floatwrite(t,nt,sout);
*/
		    /* prepare for CG */
		    fatomo_set(t,m[is]);

		    i = 0;
		    for (it=0; it < nt; it++) {
			if (m[is][it] == 1) {
			    rhs[i] = t0[is][i]-t[it];
			    i++;
			}
		    }

		    if (iter == 0)
			rhsnorm0 += cblas_snrm2(nrhs,rhs,1);
		    rhsnorm += cblas_snrm2(nrhs,rhs,1);

		    /* solve ds */
		    if (NULL == sf_getstring("topo"))
			sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt, nt,nrhs,ds,rhs,cgiter,eps,"verb",false,"end");
		    else
			sf_solver_reg(fatomo_lop,sf_cgstep,sf_igrad2_lop,2*nt, nt,nrhs,ds,rhs,cgiter,eps,"known",k,"x0",x0,"verb",false,"end");
		    
		    sf_cgstep_close();

		    /* collect gradients: modify if */
		    for (it=0; it < nt; it++) {
			if (m[is][it] != 1)
			    gs[it] += ds[it]/nshot;
		    }

		}

                rate = rhsnorm/rhsnorm0; 
		sf_warning("L2 misfit after iteration %d of %d: %g",iter,niter,rate);

		if (grad != NULL) sf_floatwrite(gs,nt,grad);
		if (norm != NULL) sf_floatwrite(&rate,1,norm);

		/* update slowness */

		for (it=0; it < nt; it++) {
		    s[it] = (s[it]+gs[it])*(s[it]+gs[it])/s[it];
     		}

		/* cheating */
/*
		for (it=0; it < nt; it++) {
		    if (gs[it] < 0.)
			s[it] = (s[it]+gs[it])*(s[it]+gs[it])/s[it];
		}
*/
	    }

	    rhsnorm = 0.;

	    /* loop over all shots */
	    for (is=0; is < nshot; is++) {
		for (it=0; it < nt; it++) {
		    t[it] = 0.;
		    ds[it] = 0.;
		}
		
		/* forward fast-marching for stencil time */
		fastmarch_init(n[2],n[1],n[0]);
		
		fastmarch(t,s,flag,plane,
			  n[2],n[1],n[0],o[2],o[1],o[0],d[2],d[1],d[0],
			  source[is][2],source[is][1],source[is][0],1,1,1,order);
		
		fastmarch_close();
		
		i = 0;
		for (it=0; it < nt; it++) {
		    if (m[is][it] == 1) {
			rhs[i] = t0[is][i]-t[it];
			i++;
		    }
		}
		
		rhsnorm += cblas_snrm2(nrhs,rhs,1);
	    }
	    
	    rate = rhsnorm/rhsnorm0;
	    sf_warning("L2 misfit after iteration %d of %d: %g",iter,niter,rate);

	    if (norm != NULL) sf_floatwrite(&rate,1,norm);

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
