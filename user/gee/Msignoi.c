/* Signal and noise separation (N-dimensional). */
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

#include "signoi.h"
#include "regrid.h"

int main(int argc, char* argv[])
{
    bool spitz, prec, verb;
    int niter, sa, na, j, dim, nx, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float *dd, *ss, eps, na0, sa0;
    char varname[6], *lagfile;
    sf_filter naa, saa;
    sf_file spef, npef, dat, signoi, slag, nlag;

    sf_init(argc,argv);
    dat = sf_input("in");
    signoi = sf_output("out");

    spef = sf_input("sfilt");
    npef = sf_input("nfilt");

    dim = sf_filedims(dat,n);
    
    if (!sf_histint(spef,"n1",&sa)) sf_error("No n1= in sfilt");
    if (!sf_histint(npef,"n1",&na)) sf_error("No n1= in nfilt");

    naa = sf_allocatehelix(na);
    saa = sf_allocatehelix(sa);

    if (NULL == (lagfile = sf_histstring(spef,"lag")) &&
	NULL == (lagfile = sf_getstring("slag"))) 
	sf_error("Need slag=");
    slag = sf_input(lagfile);
    if (!sf_histints(slag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(saa->lag,sa,slag);
    regrid(dim,m,n,saa);
    sf_fileclose(slag);

    if (NULL == (lagfile = sf_histstring(npef,"lag")) &&
	NULL == (lagfile = sf_getstring("nlag"))) 
	sf_error("Need nlag=");
    nlag = sf_input(lagfile);
    if (!sf_histints(nlag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(naa->lag,na,nlag);
    regrid(dim,m,n,naa);
    sf_fileclose(nlag);

    if (!sf_histfloat(spef,"a0",&sa0)) sa0=1.;
    if (!sf_histfloat(npef,"a0",&na0)) na0=1.;

    if (!sf_getfloat("epsilon",&eps)) sf_error("Need eps=");
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    sprintf(varname,"n%d",dim+1);
    sf_putint(signoi,varname,2);

    nx=1;
    for(j=0; j < dim; j++) {
	nx *= n[j];
    }

    dd = sf_floatalloc(nx);
    ss = sf_floatalloc(nx);

    sf_floatread(dd,nx,dat);
    sf_floatread(saa->flt,sa,spef);
    for (j=0; j < sa; j++) {
	saa->flt[j] /= sa0;
    }

    sf_floatread(naa->flt,na,npef);
    for (j=0; j < na; j++) {
	naa->flt[j] /= na0;
    }

    if (!sf_getbool("spitz",&spitz)) spitz=false;
    /* if use Spitz method */

    if (!sf_getbool("prec",&prec)) prec=false;
    /* if use preconditioning with Spitz */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    signoi_init (naa, saa, niter, nx, eps, verb);

    if (spitz) {
	if (prec) {
	    signoi1_lop  (false,false,nx,nx,dd,ss);
	} else {
	    signoi2_lop  (false,false,nx,nx,dd,ss);
	}
    } else {
	signoi_lop  (false,false,nx,nx,dd,ss);
    }

    sf_floatwrite(ss,nx,signoi);

    for (j=0; j < nx; j++) {
	dd[j] -= ss[j];
    }

    sf_floatwrite(dd,nx,signoi);


    exit (0);
}

