/* Local prediction filter (n-dimensional) with an adjoint flag. */
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

int main(int argc, char* argv[])
{
    bool verb, adj;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j, niter;
    float *d, *f, *g, mean, tmp;
    char key[6];
    sf_file dat, flt, mat;

    sf_init(argc,argv);

    dat = sf_input("basis");

    if (!sf_getbool("adj",&adj)) adj=false;
    
    if (!adj) {
	mat = sf_input("in");
	flt = sf_output("out");
    } else {
	flt = sf_input("in");
	mat = sf_output("out");
    }

    ndim = sf_filedims(dat,n);
    if(!adj){
	mdim = sf_filedims(mat,m);
	sf_putint(flt,"ndim",mdim);
	sf_putints(flt,"n",m,mdim);
    }else{
	if (!sf_histint(flt,"ndim",&mdim)) sf_error("No ndim= in input");
	if (!sf_histints(flt,"n",m,mdim)) sf_error("No n= in input");
    }

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
	if(!adj){
		snprintf(key,6,"d%d",j+1);
		sf_histfloat(dat,key,&tmp);
		sf_putfloat(flt,key,tmp);

		snprintf(key,6,"o%d",j+1);
		sf_histfloat(dat,key,&tmp);
		sf_putfloat(flt,key,tmp);
	}
	
	snprintf(key,6,"n%d",j+1);
	if(!adj){
	    sf_putint(flt,key,n[j]);
	}else{
	    sf_putint(mat,key,1);
	}
    }
    n12 = nd*ns;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    d = sf_floatalloc(n12);
    f = sf_floatalloc(n12);
    g = sf_floatalloc(nd);

    sf_multidivn_init(ns, mdim, nd, m, rect, d, NULL, verb); 
    
    sf_floatread(d,n12,dat);
    
    mean = 0.;
    for(i=0; i < n12; i++) {
	mean += d[i]*d[i];
    }
    if (mean == 0.) {
	sf_floatwrite(d,n12,flt);
	exit(0);
    }

    mean = sqrtf (mean/n12);
    
    for(i=0; i < n12; i++) {
	d[i] /= mean;
    }

    if (adj) {
	sf_floatread(f,n12,flt);
    } else {
	sf_floatread(g,nd,mat);
	for(i=0; i < nd; i++) {
	    g[i] /= mean;
	}
    }
    
    sf_multidivn_adj(!adj,g,f,niter);

    if (adj) {
	for(i=0; i < nd; i++) {
	    g[i] /= mean;
	}
	sf_floatwrite(g,nd,mat);
    } else {
	sf_floatwrite(f,n12,flt);
    }

    exit(0);
}
