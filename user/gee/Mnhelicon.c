/* Non-stationary helix convolution and deconvolution. */
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

#include "nhelix.h"
#include "nhelicon.h"
#include "npolydiv.h"
#include "regrid.h"

int main(int argc, char* argv[])
{
    int i, ia, na, nx, dim, ix, n123;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], *pch, *nh;
    float a0, *pp, *qq, *flt;
    char *lagfile;
    bool adj, div;
    nfilter aa;
    sf_file in, out, fnh, fpch, filt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    filt = sf_input("filt");

    dim = sf_filedims(in,n);

    n123 = 1;
    for (i=0; i < dim; i++) {
	n123 *= n[i];
    }
    
    if (NULL != sf_getstring("nh")) {
	fnh = sf_input("nh");
	if (!sf_histint(fnh,"n1",&nx)) sf_error("No n1= in nh");
    } else {
	fnh = NULL;
	if (!sf_histint(filt,"n2",&nx)) sf_error("No n2= in filt");
    }

    pp = sf_floatalloc(n123);
    qq = sf_floatalloc(n123);
    pch = sf_intalloc(n123);
    nh = sf_intalloc(nx);

    if (NULL != sf_getstring("pch")) {
	fpch = sf_input("pch");
	
	if (SF_INT != sf_gettype(fpch)) sf_error("Need int pch");
	sf_intread(pch,n123,fpch);
	sf_fileclose(fpch);
    } else {
	if (nx != n123) sf_error ("Wrong dimensions: %d != %d",nx,n123);
	for (ix=0; ix < nx; ix++) {
	    pch[ix]=ix;
	}
    }

    if (NULL != fnh) {
	if (SF_INT != sf_gettype(fnh)) sf_error("Need int nh");
	sf_intread(nh,nx,fnh);
    } else {
	if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
	for (ix=0; ix < nx; ix++) {
	    nh[ix] = na;
	}
    }

    aa = nallocate(nx, n123, nh, pch);
    free (pch);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;

    if (NULL == (lagfile = sf_getstring("lag")) /* file with filter lags */
	&&
	NULL == (lagfile = sf_histstring(filt,"lag"))) sf_error("Need lag=");
    lag = sf_input(lagfile);
    for (ix=0; ix < nx; ix++) {
	sf_intread(aa->hlx[ix]->lag,aa->hlx[ix]->nh,lag);
    }

    if (sf_histints (lag,"n",m,dim)) {
	for (ix=0; ix < nx; ix++) {
	    regrid (dim, m,n, aa->hlx[ix]);
	}
    }

    if (!sf_getbool ("adj",&adj)) adj=false;
    /* if y, do adjoint operation */
    if (!sf_getbool ("div",&div)) div=false;
    /* if y, do inverse operation (deconvolution) */

    if (adj) {
	sf_floatread (qq,n123,in);
    } else {
	sf_floatread (pp,n123,in);
    }

    for (ix=0; ix < nx; ix++) {
	flt = aa->hlx[ix]->flt;
	na = aa->hlx[ix]->nh;
	sf_floatread (flt,na,filt);
	for (ia=0; ia < na; ia++) {
	    flt[ia] /= a0;
	}
    }

    if (div) {
	npolydiv_init (n123, aa);
	npolydiv_lop (adj,false,n123,n123,pp,qq);
    } else {
	nhelicon_init (aa);
	nhelicon_lop (adj,false,n123,n123,pp,qq);
    }

    if (adj) {
	sf_floatwrite (pp,n123,out);
    } else {
	sf_floatwrite (qq,n123,out);
    }


    exit(0);
}

/* 	$Id: Mnhelicon.c 4796 2009-09-29 19:39:07Z ivlad $	 */
