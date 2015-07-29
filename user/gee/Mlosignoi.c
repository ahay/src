/* Local signal and noise separation (N-dimensional).

Signal and noise separation by inversion (super-deconvolution).
Uses the helix and patching technologies.
*/
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

#include "tent.h"
#include "patching.h"
#include "signoi.h"

int main(int argc, char* argv[])
{
    int n123, n1, i, ik, dim, nk, nf, sf, niter, nw;
    int n[SF_MAX_DIM], w[SF_MAX_DIM], k[SF_MAX_DIM];
    int sa[SF_MAX_DIM], na[SF_MAX_DIM], sc[SF_MAX_DIM], nc[SF_MAX_DIM];
    int ma[SF_MAX_DIM], mc[SF_MAX_DIM]; 
    float *data, *wind, *sign, eps, di, dabs;
    char varname[6], *lagfile;
    sf_filter saa, naa, sbb, nbb, scc, ncc;
    sf_file dat, signal, spef, npef, slag, nlag;

    sf_init (argc,argv);
    dat = sf_input("in");
    signal = sf_output("out");

    spef = sf_input("sfilt");
    npef = sf_input("nfilt");

    n123 = sf_filesize(dat);
    if (!sf_histint(spef,"dim",&dim)) sf_error("No dim= in sfilt");

    n1 = 1;
    for (i=0; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	if (!sf_histint(dat,varname,n+i)) 
	    sf_error("No %s= in input",varname);
	n1 *= n[i];
    }

    if (!sf_histints(spef,"w",w,dim)) sf_error("No w= in sfilt");
    if (!sf_histints(spef,"k",k,dim)) sf_error("No k= in sfilt");

    if (!sf_histints(spef,"a",sa,dim)) sf_error("No a= in sfilt");
    if (!sf_histints(npef,"a",na,dim)) sf_error("No a= in nfilt");

    if (!sf_histints(spef,"center",sc,dim)) sf_error("No center= in sfilt");
    if (!sf_histints(npef,"center",nc,dim)) sf_error("No center= in nfilt");

    nk=nw=1;
    for (i=0; i < dim; i++) {
	nw *= w[i];
	nk *= k[i];
    }

    if (!sf_histint(spef,"n1",&sf)) sf_error("No n1= in sfilt");
    if (!sf_histint(npef,"n1",&nf)) sf_error("No n1= in nfilt");

    sbb = sf_allocatehelix(sf);
    nbb = sf_allocatehelix(nf);

    if (NULL == (lagfile = sf_histstring(spef,"lag")) &&
	NULL == (lagfile = sf_getstring("slag"))) 
	sf_error("Need slag=");
    slag = sf_input(lagfile);
    if (NULL == (lagfile = sf_histstring(npef,"lag")) &&
	NULL == (lagfile = sf_getstring("nlag"))) 
	sf_error("Need nlag=");
    nlag = sf_input(lagfile);

    sf_intread(sbb->lag,sf,slag);
    sf_intread(nbb->lag,nf,nlag);

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization parameter */
    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    data = sf_floatalloc(n123);
    sign = sf_floatalloc(n123);

    sf_floatread(data,n123,dat);

    dabs = fabsf(data[0]);
    for (i=1; i < n123; i++) {
	di = fabsf(data[i]);
	if (di > dabs) dabs=di;
    }

    for (i=0; i < n123; i++) {
	data[i] /= dabs;
    }

    saa = (sf_filter) sf_alloc(nk,sizeof(*saa));
    naa = (sf_filter) sf_alloc(nk,sizeof(*naa));

    for (ik=0; ik < nk; ik++) {
	scc = saa+ik;
	ncc = naa+ik;
	scc->nh = sf;
	ncc->nh = nf;
	scc->flt = sf_floatalloc(sf);
	ncc->flt = sf_floatalloc(nf);
	scc->lag = sbb->lag;
	ncc->lag = nbb->lag;
	scc->mis = NULL;
	ncc->mis = NULL;
    }

    wind = sf_floatalloc(nw);

    for (i=0; i < dim; i++) {
	mc[i] = SF_MIN(sc[i],nc[i]);
	ma[i] = SF_MIN(sa[i],na[i]);
    }

    tent (dim, w, mc, ma, wind);
 
    for (i=0; i < n123-n1+1; i += n1) {
	signoi_init (naa, saa, niter, nw, eps, false);
	for (ik=0; ik < nk; ik++) {
	    sf_floatread((naa+ik)->flt,nf,npef);
	    sf_floatread((saa+ik)->flt,sf,spef);
	}
	patching (signoi_lop, data+i, sign+i, dim, k, n, w, wind);
    }

    sf_floatwrite (sign,n123,signal);

    exit(0);
}

/* 	$Id: Mlosignoi.c 7591 2011-08-15 17:11:55Z sfomel $	 */
