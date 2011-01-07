/* Local prediction filter for complex numbers (n-dimensional). */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "cmultidivn.h"
#include "cweight2.h"

int main(int argc, char* argv[])
{
    bool verb;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j, niter;
    sf_complex *d, *f, *g;
    float mean;
    char key[6];
    sf_file dat, flt, mat, pre;

    sf_init(argc,argv);

    dat = sf_input("in");
    mat = sf_input("match");
    flt = sf_output("out");

    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    ndim = sf_filedims(dat,n);
    mdim = sf_filedims(mat,m);

    if (mdim > ndim) 
	sf_error("Wrong dimensions: %d > %d",mdim,ndim);
 
    nd = 1;
    for (j=0; j < mdim; j++) {
	if (m[j] != n[j]) 
	    sf_error("Size mismatch [n%d]: %d != %d",
		     j+1,m[j],n[j]);
	nd *= m[j];
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1;
    }
    for (ns = 1; j < ndim; j++) {
	ns *= n[j];
	if (pre) {
	    snprintf(key,6,"n%d",j+1);
	    sf_putint(pre,key,1);
	}
    }
    n12 = nd*ns;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    d = sf_complexalloc(n12);
    f = sf_complexalloc(n12);
    g = sf_complexalloc(nd);

    cmultidivn_init(ns, mdim, nd, m, rect, d, verb); 
    
    sf_complexread(d,n12,dat);
    sf_complexread(g,nd,mat);

    mean = 0.;
    for(i=0; i < n12; i++) {
#ifdef SF_HAS_COMPLEX_H
	mean += crealf(conjf(d[i])*d[i]);
#else
	mean += crealf(sf_cmul(conjf(d[i]),d[i]));
#endif
    }
    if (mean == 0.) {
	sf_complexwrite(d,n12,flt);
	exit(0);
    }

    mean = sqrtf (n12/mean);
    
    for(i=0; i < n12; i++) {
#ifdef SF_HAS_COMPLEX_H
	d[i] *= mean;
#else
	d[i] = sf_crmul(d[i],mean);
#endif
    }
    for(i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	g[i] *= mean;
#else
	g[i] = sf_crmul(g[i],mean);
#endif
    }
    
    cmultidivn (g,f,niter);
    sf_complexwrite(f,n12,flt);

    if (pre) {
	for(i=0; i < n12; i++) {
#ifdef SF_HAS_COMPLEX_H
	    d[i] /= mean;
#else
	    d[i] = sf_crmul(d[i],1./mean);
#endif
	}
	
	cweight2_lop(false,false,n12,nd,f,g);
	sf_complexwrite(g,nd,pre);
    }


    exit(0);
}
