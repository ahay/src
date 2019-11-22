/* Streaming prediction filter */
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

#include "stmultidiv.h"

int main(int argc, char* argv[])
{
    char key[6];
    int n[SF_MAX_DIM], m[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j;
    float **d, **f, *g, mean, lam;
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
    }
    for (ns = 1; j < ndim; j++) {
	ns *= n[j];
	if (pre) {
	    snprintf(key,6,"n%d",j+1);
	    sf_putint(pre,key,1);
	}
    }
    n12 = nd*ns;

    if (!sf_getfloat("lambda",&lam)) lam=1.0f;
    /* smoothing parameter */

    d = sf_floatalloc2(nd,ns);
    f = sf_floatalloc2(nd,ns);
    g = sf_floatalloc(nd);

    stmultidiv_init(ns,nd,d,lam);
    
    sf_floatread(d[0],n12,dat);
    sf_floatread(g,nd,mat);

    mean = 0.;
    for(i=0; i < n12; i++) {
	mean += d[0][i]*d[0][i];
    }
    if (mean == 0.) {
	sf_floatwrite(d[0],n12,flt);
	exit(0);
    }

    mean = sqrtf (mean/n12);
    
    for(i=0; i < n12; i++) {
	d[0][i] /= mean;
    }
    for(i=0; i < nd; i++) {
	g[i] /= mean;
    }
    
    stmultidiv (g,f);
    sf_floatwrite(f[0],n12,flt);

    if (pre) {
	for(i=0; i < n12; i++) {
	    d[0][i] *= mean;
	}

	sf_weight2_init(ns,nd,d[0]);
	sf_weight2_lop(false,false,n12,nd,f[0],g);
	sf_floatwrite(g,nd,pre);
    }


    exit(0);
}
