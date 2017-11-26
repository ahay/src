/* Local Prediction-Error Filter (1-D, 2-D, and 3-D). */
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

#include <math.h>
#include <float.h>

#include <rsf.h>

#include "steepdip.h"
#include "createhelix.h"
#include "lopef.h"
#include "printfilter.h"
#include "bound.h"
#include "random.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], w[3], k[3], a[3], gap[3], center[3];
    int n123, n1, dim, dim1, nk, i, j, ik, na;
    bool stp;
    float *data, d[3], o[3], *mask, vel, tgap, dabs, di;
    char varname[6], *lagfile;
    sf_filter aa, bb, ak;
    sf_file dat, pef, lag, known;

    sf_init(argc,argv);
    dat = sf_input("in");
    pef = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */
    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(pef,"lag",lagfile);
    
    dim = sf_filedims(dat,n);
    if (dim > 3) sf_error("Need 3 dimensions or less");

    if (!sf_getint("dim",&dim1)) dim1=dim;
    /* PEF dimensionality */
    if (dim1 > dim) dim1=dim;
    sf_putint(pef,"dim",dim1);

    if (!sf_getints("w",w,dim1)) sf_error("Need w=");
    /* window size */
    for (j=0; j < dim1; j++) {
	if (w[j] > n[j]) w[j] = n[j];
    }
    sf_putints(pef,"w",w,dim1);

    if (!sf_getints("a",a,dim1)) sf_error("Need a=");
    /* filter size */
    sf_putints(pef,"a",a,dim1);

    if (!sf_getints("k",k,dim1)) {
	/* number of windows */
	for (j=0; j < dim1; j++) {
	    k[j] = 1.5 * n[j] / (w[j] - a[j] + 1.);
	}
    }
    sf_putints(pef,"k",k,dim1);

    if (!sf_getints("gap",gap,dim1)) {
	/* filter gap */
	for (j=0; j < dim1; j++) {
	    gap[j] = 0;
	}
    }
    sf_putints(pef,"gap",gap,dim1);

    if (!sf_getints("center",center,dim1)) {
	/* filter center */
	for (j=0; j < dim1-1; j++) {
	    center[j] = (a[j+1] > 1)? a[j]/2: 0;	    
	}
	center[dim1-1] = 0;
    }
    sf_putints(pef,"center",center,dim1);

    if (!sf_getbool("steepdip",&stp)) stp=false;
    /* if y, do steep-dip PEF estimation */

    for (j=0; j < dim1; j++) {
	sprintf(varname,"d%d",j+1);
	if (!sf_histfloat(dat,varname,d+j)) 
	    sf_error("No %s= in input",varname);
	sprintf(varname,"o%d",j+1);
	if (!sf_histfloat(dat,varname,o+j)) 
	    sf_error("No %s= in input",varname);
    }

    if (stp) {
	if (!sf_getfloat ("vel",&vel)) vel=1.7;
	/* velocity for steep-dip decon */

	if (!sf_getfloat("tgap",&tgap)) tgap=0.030;
	/* time gap for steep-dip decon */
 
	bb = steep(dim1, w, a, d, vel, tgap);
    } else {
	bb = createhelix(dim1, w, center, gap, a); 
    }

    sf_putints (lag,"n",w,dim1);

    bound (dim1, false, w, w, a, bb);
    for (i=0; i < bb->nh; i++) {
	bb->flt[i] = 2.;
    }
    print(dim1, w, center, a, bb);

    n123=n1=nk=1;
    for (j=0; j < dim; j++) {
	n123 *= n[j];
	if (j < dim1) {
	    n1 = n123;
	    nk *= k[j];
	}
    }
    na = bb->nh;

    sf_putint(pef,"n1",na);
    sf_putint(lag,"n1",na);

    for (j=0; j < dim; j++) {
	sprintf(varname,"n%d",j+2);
	sf_putint(lag,varname,1);
	if (j < dim1) {
	    sf_putint(pef,varname,k[j]);
	    sprintf(varname,"o%d",j+2);
	    sf_putfloat(pef,varname,o[j]+0.5*w[j]*d[j]);
	    sprintf(varname,"d%d",j+2);
	    sf_putfloat(pef,varname,(n[j]-w[j])/(k[j]-1.)*d[j]);
	} else if (j == dim1) {
	    sf_putint(pef,varname,n123/n1);
	} else {
	    sf_putint(pef,varname,1);
	}
    }

    sf_intwrite(bb->lag,na,lag);
    
    data = sf_floatalloc(n123);

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

    dabs = fabsf(data[0]);
    for (i=1; i < n123; i++) {
	di = fabsf(data[i]);
	if (di > dabs) dabs=di;
    }

    random_init (2004);
    for (i=0; i < n123; i++) {
	data[i] = data[i]/dabs + 100.*FLT_EPSILON*(random0()-0.5);
    }

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

    for (i=0; i < n123-n1+1; i += n1) {
	find_lopef (dim1, data+i, aa, k, n, w, mask);

	for (ik=0; ik < nk; ik++) {
	    sf_floatwrite ((aa+ik)->flt,na,pef);
	}
    }



    exit(0);
}

/* 	$Id$	 */
