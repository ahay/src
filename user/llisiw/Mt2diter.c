/* Time-to-depth conversion (linear operator) */
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

#include "t2diter.h"

int main(int argc, char* argv[])
{
    bool velocity, adj, verb, shape;
    int dim, i, j, n[3], rect[3], it, nt, nt0, nx0, cgiter;
    float d[3], o[3], dt0, dx0, ot0, ox0, eps, tol, *p=NULL, *p0=NULL;
    float *s0, *t0, *x0, *vd, *rhs, *ds, *vdt, *vdx;
    int *f0, *m0=NULL;
    char key[6], *what;
    sf_file in, out, is0=NULL, it0=NULL, ix0=NULL, if0=NULL, dix=NULL, mask=NULL, prec=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="invert";
    /* what to compute (default inversion) */

    /* read dimension */
    if (NULL == sf_getstring("s0")) sf_error("No velocity input s0=");
    is0 = sf_input("s0");
    
    dim = sf_filedims(is0,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(is0,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(is0,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; d[2] = d[1]; o[2] = o[1];
    }

    /* read velocity */
    s0 = sf_floatalloc(nt);
    sf_floatread(s0,nt,is0);
    sf_fileclose(is0);

    /* read t0 */
    if (NULL == sf_getstring("t0")) sf_error("No time input t0=");
    it0 = sf_input("t0");
    t0 = sf_floatalloc(nt);
    sf_floatread(t0,nt,it0);
    sf_fileclose(is0);

    /* read x0 */
    if (NULL == sf_getstring("x0")) sf_error("No position input x0=");
    ix0 = sf_input("x0");
    x0 = sf_floatalloc(nt);
    sf_floatread(x0,nt,ix0);
    sf_fileclose(ix0);

    /* read f0 */
    if (NULL == sf_getstring("f0")) sf_error("No flag input f0=");
    if0 = sf_input("f0");
    f0 = sf_intalloc(nt);
    sf_intread(f0,nt,if0);
    sf_fileclose(if0);

    /* read mask */
    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	m0 = sf_intalloc(nt);
	sf_intread(m0,nt,mask);
	sf_fileclose(mask);
    }

    /* read preconditioner */
    if (NULL != sf_getstring("prec")) {
	prec = sf_input("prec");
	p0 = sf_floatalloc(nt);
	sf_floatread(p0,nt,prec);
	sf_fileclose(prec);
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

    /* Dix velocity derivative in t0 */
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

    /* Dix velocity derivative in x0 */
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

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* y, inputs are velocity / n, slowness-squared */

    if (velocity) {
	for (it=0; it < nt; it++) {
	    s0[it] = 1./s0[it]*1./s0[it];
	}
    }

    /* allocate memory */
    ds = sf_floatalloc(nt);
    rhs = sf_floatalloc(nt);

    /* initialization */
    t2d_init(dim,n,d,nt0,dt0,ot0,nx0,dx0,ox0);
    
    /* set-up */
    t2d_set(t0,x0,f0,s0,vd,vdt,vdx,m0,p0);
    
    switch (what[0]) {
	case 'l': /* linear operator */

	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag */

	    /* read input */
	    if (adj)
		sf_floatread(rhs,nt,in);
	    else
		sf_floatread(ds,nt,in);
	    
	    /* apply operator */
	    if (adj)
		t2d_oper(true,false,nt,nt,ds,rhs);
	    else
		t2d_oper(false,false,nt,nt,ds,rhs);
	    
	    /* write output */
	    if (adj)
		sf_floatwrite(ds,nt,out);
	    else
		sf_floatwrite(rhs,nt,out);
	    
	    break;
	    
	case 't': /* linear prediction of dt */
	    
	    /* read input */
	    sf_floatread(ds,nt,in);
	    
	    /* compute dt */
	    t2d_dt(ds,rhs);
	    
	    /* write output */
	    sf_floatwrite(rhs,nt,out);
	    
	    break;
	    
	case 'x': /* linear prediction of dx */

	    /* read input */
	    sf_floatread(ds,nt,in);
	    
	    /* compute dt */
	    t2d_dt(ds,rhs);
	    
	    /* compute dx */
	    t2d_dx(rhs,ds);

	    /* write output */
	    sf_floatwrite(ds,nt,out);
	    
	    break;  
	    
	case 'c': /* cost */
	    
	    /* compute with upwind stencil */
	    t2d_cost(rhs);
	    
	    /* write output */
	    sf_floatwrite(rhs,nt,out);
	    
	    break;
	    
	case 'i': /* inversion */
	    
	    if (!sf_getbool("shape",&shape)) shape=false;
	    /* regularization (default Tikhnov) */
	    
	    if (!sf_getfloat("eps",&eps)) eps=0.1;
	    /* regularization parameter */

	    if (!sf_getbool("verb",&verb)) verb=false;
	    /* verbosity flag */

	    if (shape) {
		/* shaping */

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
		/* Tickhnov */

		/* initialize 2D gradient operator */
		sf_igrad2_init(n[0],n[1]);
	    }

	    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
	    /* number of CG iterations */

	    /* read input */
	    sf_floatread(rhs,nt,in);

	    /* solve ds */
	    if (shape) {
		if (p0 == NULL)
		    sf_conjgrad(NULL,t2d_oper,sf_repeat_lop,p,ds,rhs,cgiter);
		else /* shaping with data preconditioner */
		    sf_conjgrad(t2d_prec,t2d_oper,sf_repeat_lop,p,ds,rhs,cgiter);
	    } else {
		sf_solver_reg(t2d_oper,sf_cgstep,sf_igrad2_lop,2*nt,nt,nt,ds,rhs,cgiter,eps,"verb",verb,"end");
		sf_cgstep_close();
	    }

	    /* write output */
	    sf_floatwrite(ds,nt,out);

	    break;
    }
    
    exit(0);
}
