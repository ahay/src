/* First-arrival Tomography */
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

#include "fatomo.h"
#include "fastmarch.h"

int main(int argc, char* argv[])
{
    bool plane[3], **m;
    int dim, n[3], i, nm, is, ns, nrhs, iter, niter, cgiter, order, *flag, **mtemp;
    float d[3], o[3], **t, **s, *slow, *stemp, *rhs, **record;
    char key[4];
    sf_file time, mask, shot, sinp, sout, gradient;

    sf_init(argc,argv);
    sinp = sf_input("in");
    mask = sf_input("mask");
    shot = sf_input("shot");
    time = sf_input("record");
    sout = sf_output("out");
 
    if (NULL != sf_getstring("gradient")) {
	gradient = sf_output("gradient");
    } else {
	gradient = NULL;
    }

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    if (!sf_getint("cgiter",&cgiter)) cgiter=200;
    /* number of conjugate gradient iterations */
    
    if (!sf_getint("order",&order)) order=2;
    /* fast marching accuracy order */

    /* scan for model dimensions */
    dim = sf_filedims(sinp,n);

    nm = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(sinp,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(sinp,key,o+i)) o[i]=0.;
	nm *= n[i];
	plane[i] = false;
    }

    /* 2D case */
    if (dim < 3) {
	n[2] = 1;
	o[2] = 0.;
	d[2] = d[1];
	plane[2] = false;
    }

    /* scan for shot dimensions */
    if (!sf_histint(shot,"n2",&ns)) ns=1;

    s     = sf_floatalloc2(3,ns);
    t     = sf_floatalloc2(nm,ns);
    m     = sf_boolalloc2(nm,ns);
    slow  = sf_floatalloc(nm);
    flag  = sf_intalloc(nm);
    mtemp = sf_intalloc2(nm,ns);
    stemp = sf_floatalloc(nm);

    if (NULL != gradient) {
	if (dim < 3) sf_putint(gradient,"n3",niter);
	else         sf_putint(gradient,"n4",niter);
    }

    sf_floatread(s[0],3*ns,shot);
    sf_intread(mtemp[0],nm*ns,mask);
    sf_floatread(slow,nm,sinp);

    nrhs = 0;
    for (is=0; is < ns; is++) {
	for (i=0; i < nm; i++) {
	    if (mtemp[is][i] == 1) {
		m[is][i] = true;
		nrhs++;
	    } else {
		m[is][i] = false;
	    }
	}
    }

    rhs = sf_floatalloc(nrhs);
    record = sf_floatalloc2(nrhs/ns,ns);
    sf_floatread(record[0],nrhs,time);

    fatomo_init(ns,dim,n,d,m);
    fastmarch_init(n,o,d,order);

    for (iter=0; iter < niter; iter++) {

	sf_warning("iteration %d of %d",iter+1,niter);

	/* fast marching forward eikonal */    
	for (is=0; is < ns; is++) {
	    fastmarch(t[is],slow,s[is]);
	}
	
	/* set stencil and right-hand side */
	fatomo_set(t,record,rhs);

	sf_solver(fatomo_lop,sf_cgstep,nm,nrhs,stemp,rhs,cgiter,"verb",false,"end");
	
	/* update slowness */
	for (i=0; i < nm; i++) {
	    slow[i] = slow[i]+2*stemp[i];
	    if (slow[i] <= 0.) slow[i] = FLT_EPSILON;
	}

	if (NULL != gradient) sf_floatwrite(stemp,nm,gradient);
    }

    sf_floatwrite(slow,nm,sout);

    exit(0);
}
