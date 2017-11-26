/* Local Matched filter for coherent noise removal (1-D, 2-D, and 3-D). */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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

#include <math.h>
#include <float.h>

#include <rsf.h>
#include <rsfgee.h>

#include "lomatch.h"

int main(int argc, char* argv[])
{
    bool verb;
	int n[SF_MAX_DIM], w[3], k[3], a[3], gap[3], center[3];
    int n123, n1, dim, dim1, nk, i, j, ik, na;
    int niter,liter; float perc;
    float *data=NULL, *match=NULL, d[3], o[3], *mask;
    //float dabs, di, mabs, mi;
    float eps;
    char varname[6], *lagfile;
    sf_filter aa, bb, ak;
    sf_file dat=NULL, mat=NULL, mcf=NULL, lag=NULL, known=NULL;

    sf_init(argc,argv);
    dat = sf_input("in");
    mcf = sf_output("out");

    if (NULL != sf_getstring("match")) mat=sf_input("match");
    	else sf_error("No match file specified");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */
    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(mcf,"lag",lagfile);
    
    dim = sf_filedims(dat,n);
    if (dim > 3) sf_error("Need 3 dimensions or less");

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of POCS iterations: =1 L2 norm ; >1 L1 norm */

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening [L1 norm]*/

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* dumping parameter x=(M'M+eps*I)M'd */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("dim",&dim1)) dim1=dim;
    /* matched filter  dimensionality */
    if (dim1 > dim) dim1=dim;
    sf_putint(mcf,"dim",dim1);

    if (!sf_getints("w",w,dim1)) sf_error("Need w=");
    /* window size */
    for (j=0; j < dim1; j++) {
    	if (w[j] > n[j]) w[j] = n[j];
    }
    sf_putints(mcf,"w",w,dim1);

    if (!sf_getints("a",a,dim1)) sf_error("Need a=");
    /* filter size */
    sf_putints(mcf,"a",a,dim1);

    if (!sf_getints("k",k,dim1)) {
  	/* number of windows */
    	for (j=0; j < dim1; j++) {
    		k[j] = 1.5 * n[j] / (w[j] - a[j] + 1.);
		}
    }
    sf_putints(mcf,"k",k,dim1);

    if (!sf_getints("gap",gap,dim1)) {
  	/* filter gap */
  	for (j=0; j < dim1; j++) {
  	    gap[j] = 0;
  	}
    }
    sf_putints(mcf,"gap",gap,dim1);

    if (!sf_getints("center",center,dim1)) {
  	/* filter center */
  	for (j=0; j < dim1-1; j++) {
  	    center[j] = (a[j+1] > 1)? a[j]/2: 0;
  	}
  	center[dim1-1] = 0;
    }
    sf_putints(mcf,"center",center,dim1);

    for (j=0; j < dim1; j++) {
	sprintf(varname,"d%d",j+1);
	if (!sf_histfloat(dat,varname,d+j)) 
	    sf_error("No %s= in input",varname);
	sprintf(varname,"o%d",j+1);
	if (!sf_histfloat(dat,varname,o+j)) 
	    sf_error("No %s= in input",varname);
    }


    bb = createhelix(dim1, w, center, gap, a);



    sf_putints (lag,"n",w,dim1);


    bound (dim1, false, w, w, a, bb);
    for (i=0; i < bb->nh; i++) {
	bb->flt[i] = 2.;
    }
    //print(dim1, w, center, a, bb);

    n123=n1=nk=1;
    for (j=0; j < dim; j++) {
	n123 *= n[j];
	if (j < dim1) {
	    n1 = n123;
	    nk *= k[j];
	}
    }
    na = bb->nh;

    sf_putint(mcf,"n1",na);
    sf_putint(lag,"n1",na);

    for (j=0; j < dim; j++) {
	sprintf(varname,"n%d",j+2);
	sf_putint(lag,varname,1);
	if (j < dim1) {
	    sf_putint(mcf,varname,k[j]);
	    sprintf(varname,"o%d",j+2);
	    sf_putfloat(mcf,varname,o[j]+0.5*w[j]*d[j]);
	    sprintf(varname,"d%d",j+2);
	    sf_putfloat(mcf,varname,(n[j]-w[j])/(k[j]-1.)*d[j]);
	} else if (j == dim1) {
	    sf_putint(mcf,varname,n123/n1);
	} else {
	    sf_putint(mcf,varname,1);
	}
    }

    sf_intwrite(bb->lag,na,lag);
    
    data = sf_floatalloc(n123);
    match = sf_floatalloc(n123);

    aa = (sf_filter) sf_alloc(nk,sizeof(*aa));

    if (NULL != sf_getstring ("mask")) {
	known = sf_input("mask");
	mask = sf_floatalloc(n123);
	sf_floatread(mask,n123,known);
	sf_fileclose(known);
    } else {
	mask = NULL;
    }

    sf_floatread(data,n123,dat);
    sf_floatread(match,n123,mat);


    for (ik=0; ik < nk; ik++) {
	ak = aa+ik;
	ak->nh = na;
	ak->flt = sf_floatalloc(na);
	for (i=0; i < na; i++) {
	    ak->flt[i] = 0.;
	}
	ak->lag = bb->lag;
	ak->mis = bb->mis;
    }

    if (!sf_getint("liter",&liter)) liter=aa->nh;
    /* number of linear iteration*/

    for (i=0; i < n123-n1+1; i += n1) {
	find_lomatch (dim1, data+i, match+i, aa, k, n, w, mask, eps, liter,niter, perc, verb);

	for (ik=0; ik < nk; ik++) {
	    sf_floatwrite ((aa+ik)->flt,na,mcf);
	}
    }

    sf_deallocatehelix(bb);

    exit(0);
}

