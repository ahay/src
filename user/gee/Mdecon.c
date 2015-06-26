/* Deconvolution (N-dimensional).

Uses the helix and patching technology.
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

#include <math.h>

#include <rsf.h>

#include "tent.h"
#include "patching.h"
#include "loconvol.h"

int main(int argc, char* argv[])
{
    int n123, n1,i,ik, dim,nk,nw, na, r1, n[SF_MAX_DIM], w[SF_MAX_DIM];
    int k[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM];
    float *data, *wind, *resi, *sign, signvar, datavar, di, dabs;
    char varname[6], *lagfile;
    bool predictive;
    sf_filter aa, bb, cc;
    sf_triangle tr=NULL;
    sf_file in, out, filt, lag;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    filt = sf_input("filt");

    n123 = sf_filesize(in);
    if (!sf_histint(filt,"dim",&dim)) sf_error("No dim= in filt");

    n1 = 1;
    for (i=0; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	if (!sf_histint(in,varname,n+i)) 
	    sf_error("No %s= in input",varname);
	n1 *= n[i];
    }

    if (!sf_histints(filt,"w",w,dim)) sf_error("No w= in sfilt");
    if (!sf_histints(filt,"k",k,dim)) sf_error("No k= in sfilt");
    if (!sf_histints(filt,"a",a,dim)) sf_error("No a= in sfilt");
    if (!sf_histints(filt,"center",center,dim)) 
	sf_error("No center= in sfilt");

    nk=nw=1;
    for (i=0; i < dim; i++) {
	nw *= w[i];
	nk *= k[i];
    }

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in sfilt");
    bb = sf_allocatehelix(na);
    
    if (NULL == (lagfile = sf_histstring(filt,"lag")) &&
	NULL == (lagfile = sf_getstring("lag"))) 
	sf_error("Need lag=");
    lag = sf_input(lagfile);
    sf_intread(bb->lag,na,lag);

    if (!sf_getbool ("predictive",&predictive)) predictive=false;
    /* if y, do predictive deconvolution */

    if (!sf_getint("rect1",&r1)) r1=0;
    /* smoothing in the first axis */

    data = sf_floatalloc(n123);
    resi = sf_floatalloc(n123);
    sign = sf_floatalloc(n123);

    sf_floatread (data,n123,in);

    dabs = fabsf(data[0]);
    for (i=1; i < n123; i++) {
	di = fabsf(data[i]);
	if (di > dabs) dabs=di;
    }

    for (i=0; i < n123; i++) {
	data[i] /= dabs;
    }

    aa = (sf_filter) sf_alloc(nk,sizeof(*aa));

    for (ik=0; ik < nk; ik++) {
	cc = aa+ik;
	cc->nh = na;
	cc->flt = sf_floatalloc(na);
	cc->lag = bb->lag;
	cc->mis = NULL;
    }

    wind = sf_floatalloc(nw);
    tent (dim, w, center, a, wind);

    for (i=0; i < n123-n1+1; i += n1) {
	loconvol_init (aa);
	for (ik=0; ik < nk; ik++) {     
	    sf_floatread ((aa+ik)->flt,na,filt);
	}
	patching (loconvol_lop, data+i, resi+i, dim, k, n, w, wind);
    }

    if (predictive) {
	signvar = 0.;
	datavar = 0.;
	for (i=0; i < n123; i++) {
	    sign[i] = data[i] - resi[i];
	    signvar += sign[i]*sign[i];
	    datavar += data[i]*data[i];
	}
	signvar = sqrtf(signvar/datavar);
	for (i=0; i < n123; i++) {
	    if (resi[i] == 0.) sign[i] = data[i]*signvar;
	}
	sf_floatwrite(sign,n123,out);
    } else if (r1 > 0) {
	tr = sf_triangle_init(r1,n[0],false);
	for (i=0; i < n123-n[0]+1; i += n[0]) {
	    sf_smooth2(tr,0,1,false,resi+i);
	}
	sf_floatwrite(resi,n123,out);
    } else {
	signvar = 0.;
	datavar = 0.;
	for (i=0; i < n123; i++) {
	    signvar += resi[i]*resi[i];
	    datavar += data[i]*data[i];
	}
	signvar = sqrtf(signvar/datavar);
	for (i=0; i < n123; i++) {
	    if (resi[i] == 0.) resi[i] = data[i]*signvar;
	}
	sf_floatwrite(resi,n123,out);
    }

    exit(0);
}

/* 	$Id$	 */
