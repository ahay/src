/* Iterative time-to-depth velocity conversion */
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

#include "irays.h"
#include "t2diter.h"

int main(int argc, char* argv[])
{
    bool velocity, verb, shape;
    int dim, i, j, n[3], rect[3], it, nt, order, nt0, nx0;
    int iter, niter, iline, nline, cgiter, *f0, *m0=NULL;
    float d[3], o[3], dt0, dx0, ot0, ox0, eps, tol, *p=NULL, *p0=NULL, thres;
    float *vd, *vdt, *vdx, *s, *t0, *x0, *ds, *rhs, *rhs0, *rhs1=NULL, error0, error1, error, scale;
    char key[6];
    sf_file in, out, dix, t_0=NULL, x_0=NULL, f_0=NULL, grad=NULL, cost=NULL, mini=NULL, prec=NULL;

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
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* y, input is velocity / n, slowness-squared */

    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./s[it]*1./s[it];
	}
    }

    /* read Dix velocity */
    if (NULL == sf_getstring("dix")) sf_error("No Dix input dix=");
    dix = sf_input("dix");

    if(!sf_histint(dix,"n1",&nt0)) sf_error("No n1= in dix");
    if(!sf_histint(dix,"n2",&nx0)) sf_error("No n2= in dix");
    if(!sf_histfloat(dix,"d1",&dt0)) sf_error("No d1= in dix");
    if(!sf_histfloat(dix,"d2",&dx0)) sf_error("No d2= in dix");
    if(!sf_histfloat(dix,"o1",&ot0)) sf_error("No o1= in dix");
    if(!sf_histfloat(dix,"o2",&ox0)) sf_error("No o2= in dix");

    vd = sf_floatalloc(nt0*nx0);
    sf_floatread(vd,nt0*nx0,dix);
    sf_fileclose(dix);

    /* Dix velocity derivative in t0 (2nd order FD) */
    vdt = sf_floatalloc(nt0*nx0);
    for (i=0; i < nt0; i++) {
	for (j=0; j < nx0; j++) {
	    if (i == 0)
		vdt[j*nt0+i] = (-vd[j*nt0+i+2]+4.*vd[j*nt0+i+1]-3.*vd[j*nt0+i])/(2.*dt0);
	    else if (i == nt0-1)
		vdt[j*nt0+i] = (3.*vd[j*nt0+i]-4.*vd[j*nt0+i-1]+vd[j*nt0+i-2])/(2.*dt0);
	    else
		vdt[j*nt0+i] = (vd[j*nt0+i+1]-vd[j*nt0+i-1])/(2.*dt0);
	}
    }

    /* Dix velocity derivative in x0 (2nd order FD) */
    vdx = sf_floatalloc(nt0*nx0);
    for (j=0; j < nx0; j++) {
	for (i=0; i < nt0; i++) {
	    if (j == 0)
		vdx[j*nt0+i] = (-vd[(j+2)*nt0+i]+4.*vd[(j+1)*nt0+i]-3.*vd[j*nt0+i])/(2.*dx0);
	    else if (j == nx0-1)
		vdx[j*nt0+i] = (3.*vd[j*nt0+i]-4.*vd[(j-1)*nt0+i]+vd[(j-2)*nt0+i])/(2.*dx0);
	    else
		vdx[j*nt0+i] = (vd[(j+1)*nt0+i]-vd[(j-1)*nt0+i])/(2.*dx0);
	}
    }

    if (!sf_getint("order",&order)) order=1;
    /* fastmarch accuracy order */

    if (!sf_getfloat("thres",&thres)) thres=10.;
    /* thresholding for caustics */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of nonlinear updates */

    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
    /* number of CG iterations */

    if (!sf_getbool("shape",&shape)) shape=false;
    /* regularization (default Tikhnov) */

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* regularization parameter */

    if (!sf_getint("nline",&nline)) nline=0;
    /* maximum number of line search (default turned-off) */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

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
	
	sf_conjgrad_init(nt,nt,nt,nt,eps,tol,verb,false);
	p = sf_floatalloc(nt);
    } else {
	/* initialize 2D gradient operator */
	sf_igrad2_init(n[0],n[1]);
    }
    
    /* allocate memory for fastmarch */
    t0 = sf_floatalloc(nt);
    x0 = sf_floatalloc(nt);
    f0 = sf_intalloc(nt);

    /* allocate memory for update */
    ds  = sf_floatalloc(nt);
    rhs = sf_floatalloc(nt);

    /* output transformation matrix */
    if (NULL != sf_getstring("t0")) {
	t_0 = sf_output("t0");
	sf_putint(t_0,"n3",niter+1);
    }
    if (NULL != sf_getstring("x0")) {
	x_0 = sf_output("x0");
	sf_putint(x_0,"n3",niter+1);
    }

    /* output auxiliary label */
    if (NULL != sf_getstring("f0")) {
	f_0 = sf_output("f0");
	sf_settype(f_0,SF_INT);
	sf_putint(f_0,"n3",niter+1);
    }

    /* output gradient */
    if (NULL != sf_getstring("grad")) {
	grad = sf_output("grad");
	sf_putint(grad,"n3",niter);
    }

    /* output cost */
    if (NULL != sf_getstring("cost")) {
	cost = sf_output("cost");
	sf_putint(cost,"n3",niter+1);
    }

    /* read mask (desired minimum) */
    m0 = sf_intalloc(nt);

    if (NULL != sf_getstring("mask")) {
	mini = sf_input("mask");
	sf_intread(m0,nt,mini);
	sf_fileclose(mini);
    } else {
	for (it=0; it < nt; it++) m0[it] = -1;
    }

    /* read cost (desired minimum) */
    rhs0 = sf_floatalloc(nt);

    if (NULL != sf_getstring("mval")) {
	mini = sf_input("mval");
	sf_floatread(rhs0,nt,mini);
	sf_fileclose(mini);
    } else {
	for (it=0; it < nt; it++) rhs0[it] = 0.;
    }

    /* read preconditioner */
    if (NULL != sf_getstring("prec")) {
	prec = sf_input("prec");
	p0 = sf_floatalloc(nt);
	sf_floatread(p0,nt,prec);
	sf_fileclose(prec);

	rhs1 = sf_floatalloc(nt);
    }

    /* fastmarch initialization */
    fastmarch_init(n,o,d,order);

    /* update initialization */
    t2d_init(dim,n,d,nt0,dt0,ot0,nx0,dx0,ox0);

    /* fastmarch */
    fastmarch(t0,x0,f0,s);

    /* caustic region (2D) */
    t2d_caustic(x0,f0,n,d,thres);

    /* set up operator */
    t2d_set(t0,x0,f0,s,vd,vdt,vdx,m0,p0);

    /* evaluate cost */
    t2d_cost(rhs);

    for (it=0; it < nt; it++) {
	if (f0[it] >= 0 || m0[it] >= 0)
	    rhs[it] = 0.;
	else
	    rhs[it] -= rhs0[it];
    }

    if (p0 == NULL) {
	error0 = error1 = cblas_snrm2(nt,rhs,1);
    } else {
	for (it=0; it < nt; it++) rhs1[it] = p0[it]*rhs[it];
	error0 = error1 = cblas_snrm2(nt,rhs1,1);
    }

    /* write optional outputs */    
    if (NULL!=t_0)  sf_floatwrite(t0,nt,t_0);
    if (NULL!=x_0)  sf_floatwrite(x0,nt,x_0);
    if (NULL!=f_0)  sf_intwrite(f0,nt,f_0);
    if (NULL!=cost) sf_floatwrite(rhs,nt,cost);

    sf_warning("Start conversion, cost %g",1.);

    /* nonlinear loop */
    for (iter=0; iter < niter; iter++) {
	
	/* solve ds */
	if (shape) {
	    if (p0 == NULL)
		sf_conjgrad(NULL,t2d_oper,sf_repeat_lop,p,ds,rhs,cgiter);
	    else
		sf_conjgrad(t2d_prec,t2d_oper,sf_repeat_lop,p,ds,rhs,cgiter);
	} else {
	    sf_solver_reg(t2d_oper,sf_cgstep,sf_igrad2_lop,2*nt,nt,nt,ds,rhs,cgiter,eps,"verb",verb,"end");
	    sf_cgstep_close();
	}

	/* add ds */
	for (it=0; it < nt; it++) {
	    s[it] = s[it]+ds[it]+0.25*ds[it]*ds[it]/s[it];
	}

	/* fastmarch */
	fastmarch(t0,x0,f0,s);

	/* caustic region (2D) */
	t2d_caustic(x0,f0,n,d,thres);

	/* set up operator */
	t2d_set(t0,x0,f0,s,vd,vdt,vdx,m0,p0);

	/* evaluate cost */
	t2d_cost(rhs);

	for (it=0; it < nt; it++) {
	    if (f0[it] >= 0 || m0[it] >= 0)
		rhs[it] = 0.;
	    else
		rhs[it] -= rhs0[it];
	}
	
	if (p0 == NULL) {
	    error = cblas_snrm2(nt,rhs,1);
	} else {
	    for (it=0; it < nt; it++) rhs1[it] = p0[it]*rhs[it];
	    error = cblas_snrm2(nt,rhs1,1);
	}

	error = cblas_snrm2(nt,rhs,1);

	/* line search */
	if (nline > 0 && error >= error1) {

	    scale = 0.5;
	    for (iline=0; iline < nline; iline++) {

		for (it=0; it < nt; it++) {
		    s[it] = s[it]+(scale*ds[it])+0.25*(scale*ds[it])*(scale*ds[it])/s[it];
		}
		
		fastmarch(t0,x0,f0,s);
		t2d_caustic(x0,f0,n,d,thres);

		t2d_set(t0,x0,f0,s,vd,vdt,vdx,m0,p0);
		t2d_cost(rhs);

		for (it=0; it < nt; it++) {
		    if (f0[it] >= 0 || m0[it] >= 0)
			rhs[it] = 0.;
		    else
			rhs[it] -= rhs0[it];
		}
		
		if (p0 == NULL) {
		    error = cblas_snrm2(nt,rhs,1);
		} else {
		    for (it=0; it < nt; it++) rhs1[it] = p0[it]*rhs[it];
		    error = cblas_snrm2(nt,rhs1,1);
		}
		
		error = cblas_snrm2(nt,rhs,1);
		if (error < error1) {
		    sf_warning("Exist line search %d of %d",iline+1,nline);
		} else {
		    scale *= 0.5;
		}
	    }
	}

	error1 = error;

	/* write optional outputs */
	if (NULL!=t_0)  sf_floatwrite(t0,nt,t_0);
	if (NULL!=x_0)  sf_floatwrite(x0,nt,x_0);
	if (NULL!=f_0)  sf_intwrite(f0,nt,f_0);
	if (NULL!=cost) sf_floatwrite(rhs,nt,cost);
	if (NULL!=grad) sf_floatwrite(ds,nt,grad);

	sf_warning("Cost after iteration %d: %g",iter+1,error/error0);
    }

    /* write output */
    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./sqrtf(s[it]);
	}
    }

    sf_floatwrite(s,nt,out);

    exit(0);
}
