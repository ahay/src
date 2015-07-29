/* Multiscale missing data interpolation (N-dimensional). */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "mshelix.h"
#include "msmis2.h"
#include "bound.h"

int main(int argc, char* argv[])
{
    int na,ia, niter, id, is, ns, dim, nd, *jump;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], a[SF_MAX_DIM], b[SF_MAX_DIM];
    float a0, *mm, *kk;
    bool *known;
    msfilter msaa;
    sf_filter aa;
    char *lagfile;
    sf_file in, out, filt, lag, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    dim = sf_filedims(in,n);
    sf_putint (out,"dim",dim);

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */
    
    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");

    if (NULL == (lagfile = sf_histstring(filt,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
	/* optional input file with filter lags */
	lagfile = "lag";
    }

    lag = sf_input(lagfile);
    if (SF_INT != sf_gettype(lag)) 
	sf_error("Need int data in lag file '%s'",lagfile);
    if (!sf_histints(lag,"n",m,dim)) sf_error("No n= in lag");

    if (!sf_histint(lag,"n2",&ns)) sf_error("No n2= in lag");

    msaa = msallocate (na, ns);
    sf_intread(msaa->lag[0],na*ns,lag);
 
    jump = sf_intalloc(ns);
    if (!sf_histints(lag,"jump",jump,ns)) sf_error("No jump= in lag");

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;

    if (!sf_histints(filt,"a",a,dim)) {
	for (id=0; id < dim; id++) {
	    a[id]=1;
	}
    }

    nd=1;
    for (id=0; id < dim; id++) {
	nd *= n[id];
    }

    msaa->mis = sf_boolalloc2(nd,ns);

    aa = sf_allocatehelix(na);

    for (is=0; is < ns; is++) {
	for (ia=0; ia < na; ia++) {
	    aa->lag[ia] = msaa->lag[0][ia]*jump[is]/jump[0];
	}
	for (id=0; id < dim; id++) {
	    b[id] = a[id]*jump[is];
	}
	bound (dim, m, n, b, aa);
	for (id=0; id < nd; id++) {
	    msaa->mis[is][id] = aa->mis[id];
	}
	free (aa->mis);
	aa->mis = NULL;
    }

    sf_deallocatehelix(aa);

    mm = sf_floatalloc(nd);
    kk = sf_floatalloc(nd);
    known = sf_boolalloc(nd);

    sf_floatread(mm,nd,in);

    if (NULL != sf_getstring("mask")) {
	/* optional input mask file for known data */
	mask = sf_input("mask");
	sf_floatread(kk,nd,mask);
	sf_fileclose(mask);
	
	for (id=0; id < nd; id++) {
	    known[id] = (bool) (kk[id] != 0.);
	}
    } else {
	for (id=0; id < nd; id++) {
	    known[id] = (bool) (mm[id] != 0.);
	}
    }

    sf_floatread(msaa->flt,na,filt);
    sf_fileclose(filt);

    for (ia=0; ia < na; ia++) {
	msaa->flt[ia] /= a0;
    }

    msmis2 (niter, nd, ns, mm, msaa, known);

    sf_floatwrite(mm,nd,out);


    exit (0);
}

/* 	$Id: Mmsmiss.c 4796 2009-09-29 19:39:07Z ivlad $	 */
