/* Linearized Complex Eikonal Equation Solver */
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
#include "cpxeikonal.h"
#include "fastmarch.h"

int main(int argc, char* argv[])
{
    bool adj, velocity, plane[3], *m;
    int n[SF_MAX_DIM], it, nt, i, dim, order, *flag, cgiter, iter, niter;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *real, *imag, *t, *s, *temp, *w, *dimag, *rhs, *ref, *x0, source[3];
    float rhsnorm0, rhsnorm, rate;
    char key[4], *what;
    upgrad upgreal, upgimag;
    sf_file in, out, realin, imagin, imagref, realout, imagout, norm;

    sf_init(argc,argv);
    in      = sf_input("in");
    realin  = sf_input("realin");
    imagin  = sf_input("imagin");
    out     = sf_output("out");

    dim = sf_filedims(realin,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(realin,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(realin,key,o+i)) o[i]=0.;
	nt *= n[i]; plane[i] = false;
    }

    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1]; plane[2] = false;
    }

    if (NULL == (what = sf_getstring("what"))) what="solve";
    /* what to compute (default tomography) */

    s    = sf_floatalloc(nt);
    t    = sf_floatalloc(nt);
    real = sf_floatalloc(nt);
    imag = sf_floatalloc(nt);
    temp = sf_floatalloc(nt);

    sf_floatread(s,nt,in);

    sf_floatread(real,nt,realin);
    sf_fileclose(realin);

    sf_floatread(imag,nt,imagin);
    sf_fileclose(imagin);

    switch (what[0]) {
	case 'l': /* linear operator */

	    realout = NULL;
	    imagout = NULL;
	    imagref = NULL;

	    upgreal = upgrad_init(dim,n,d);
	    upgimag = upgrad_init(dim,n,d);
	    
	    upgrad_set(upgreal,real);
	    upgrad_set(upgimag,imag);
	    
	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */
	    
	    if (adj) {
		upgrad_adj(upgimag,temp,s);
		upgrad_inverse(upgreal,t,temp,NULL);
		upgrad_adj(upgimag,temp,t);
		
		upgrad_adj(upgreal,t,s);
		
		for (it=0; it < nt; it++) {
		    t[it] += temp[it];
		}
	    } else {
		upgrad_forw(upgimag,s,temp);
		upgrad_solve(upgreal,temp,t,NULL);
		upgrad_forw(upgimag,t,temp);
		
		upgrad_forw(upgreal,s,t);
		
		for (it=0; it < nt; it++) {
		    t[it] += temp[it];
		}
	    }
	    
	    sf_floatwrite(t,nt,out);
	    
	    break;

	case 's': /* solve complex eikonal */
	    
	    upgreal = NULL;
	    upgimag = NULL;

	    imagref = sf_input("imagref");
	    realout = sf_output("realout");
	    imagout = sf_output("imagout");

	    w     = sf_floatalloc(nt);
	    m     = sf_boolalloc(nt);
	    x0    = sf_floatalloc(nt);
	    ref   = sf_floatalloc(nt);
	    rhs   = sf_floatalloc(nt);
	    flag  = sf_intalloc(nt);
	    dimag = sf_floatalloc(nt);

	    /* read reference imag time */
	    sf_floatread(ref,nt,imagref);
	    sf_fileclose(imagref);

	    if (!sf_getint("order",&order)) order=2;
	    /* fast marching accuracy order */

	    if (!sf_getbool("velocity",&velocity)) velocity=true;
	    /* if y, the input is velocity; n, slowness squared */

	    if (!sf_getint("niter",&niter)) niter=1;
	    /* number of iterations */

	    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
	    /* number of conjugate gradient iterations */

	    if (!sf_getfloat("xshot",source+2)) source[2]=o[2]+0.5*(n[2]-1)*d[2];
	    if (!sf_getfloat("yshot",source+1)) source[1]=o[1]+0.5*(n[1]-1)*d[1];
	    if (!sf_getfloat("zshot",source))   source[0]=0.;
	    /* shot location */

	    /* output L2 norm at each iteration */
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
	    
	    /* output gradient at each iteration */
	    sf_putint(out,"n3",n[2]);
	    sf_putfloat(out,"d3",d[2]);
	    sf_putfloat(out,"o3",o[2]);
	    sf_putint(out,"n4",niter);

	    if (velocity) {
		for (it=0; it < nt; it++) {
		    s[it] = 1./s[it]*1./s[it];
		}
	    }

	    for (it=0; it < nt; it++) {
		m[it] = false;
	    }
	    for (it=0; it < n[1]; it++) {
		m[it*n[0]] = true;
	    }

	    /* initialize stencils */
	    cpxeikonal_init(dim,n,d,nt);

	    /* compute pseudo slowness */
	    cpxeikonal_set(false,imag);
/*
	    cpxeikonal_pseudo(imag,w);

	    for (it=0; it < nt; it++) {
		temp[it] = s[it]+w[it];
	    }

	    fastmarch_init(n[2],n[1],n[0]);
	    
	    fastmarch(real,temp,flag,plane,
		      n[2],n[1],n[0],
		      o[2],o[1],o[0],
		      d[2],d[1],d[0],
		      source[2],source[1],source[0],
		      1,1,1,order);
	    
	    fastmarch_close();
*/
	    /* update real stencil */
	    cpxeikonal_set(true,real);

	    /* compute right-hand side */
	    cpxeikonal_rhs(real,imag,rhs);	    
	    for (it=0; it < nt; it++) {
		rhs[it] = -rhs[it]/2;
	    }

	    rhsnorm0 = cblas_snrm2(nt,rhs,1);
	    rhsnorm = rhsnorm0;
	    rate = rhsnorm/rhsnorm0;
	    if (norm != NULL) sf_floatwrite(&rate,1,norm);

	    for (iter=0; iter < niter; iter++) {
		sf_warning("right-hand side L2 norm = %g at iter = %d",rate,iter);
		
		/* compute reference */
		for (it=0; it < nt; it++) {
		    x0[it] = ref[it]-imag[it];
		}
		
		/* solve for gradient */
		sf_solver(cpxeikonal_loop,sf_cgstep,nt,nt,dimag,rhs,cgiter,"known",m,"x0",x0,"verb",false,"end");
		sf_cgstep_close();
		
		/* update real and imag times */
		for (it=0; it < nt; it++) {
		    imag[it] += dimag[it];
/*
		    if (imag[it] <= 0.) imag[it] = FLT_EPSILON;
		    if (it < n[0]) imag[it] = 0.;
*/
		}

		cpxeikonal_set(false,imag);
		cpxeikonal_pseudo(imag,w);
		
		for (it=0; it < nt; it++) {
		    temp[it] = s[it]+w[it];
		}
		
		fastmarch_init(n[2],n[1],n[0]);
		
		fastmarch(real,temp,flag,plane,
			  n[2],n[1],n[0],
			  o[2],o[1],o[0],
			  d[2],d[1],d[0],
			  o[2]+source[2]*d[2],
			  o[1]+source[1]*d[1],
			  o[0]+source[0]*d[0],
			  1,1,1,order);
		
		fastmarch_close();

		cpxeikonal_set(true,real);
		
		cpxeikonal_rhs(real,imag,rhs);
		for (it=0; it < nt; it++) {
		    rhs[it] = -rhs[it]/2;
		}

		rhsnorm = cblas_snrm2(nt,rhs,1);
		rate = rhsnorm/rhsnorm0;
		
		if (norm != NULL) sf_floatwrite(&rate,1,norm);
		sf_floatwrite(dimag,nt,out);
	    }

	    sf_warning("right-hand side L2 norm = %g at iter = %d",rate,iter);

	    sf_floatwrite(real,nt,realout);
	    sf_floatwrite(imag,nt,imagout);

	    break;
    }

    exit(0);
}
