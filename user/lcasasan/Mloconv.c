/* Local internal and transient convolution (N-dimensional).
	Uses the helix and patching technology.
*/

/*
Copyright (C) 2010 Politecnico

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
#include <rsfgee.h>

#include "loconvol_internal.h"
#include "loconvol_transient.h"

int main(int argc, char* argv[])
{
    int n123, n1,i,ik, dim,nk,nw, na, n[SF_MAX_DIM], w[SF_MAX_DIM];
    int k[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM];
    float *data, *wind, *resi;
    char varname[6], *lagfile;
    sf_filter aa, bb, cc;
    sf_file in, out, filt, lag;
    bool transient;

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

    if (!sf_getbool("transient",&transient)) transient = false;
     /* transient convolution y/n */

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


    data = sf_floatalloc(n123);
    resi = sf_floatalloc(n123);


    sf_floatread (data,n123,in);



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
		transient ? loconvol_transient_init (aa,nw) : loconvol_internal_init (aa);
		for (ik=0; ik < nk; ik++) {
			sf_floatread ((aa+ik)->flt,na,filt);
		}
		if (transient) {
			patching (loconvol_transient_lop, data+i, resi+i, dim, k, n, w, wind);
			loconvol_transient_close();
		} else {
			patching (loconvol_internal_lop , data+i, resi+i, dim, k, n, w, wind);
		}
    }
	sf_floatwrite(resi,n123,out);


    exit(0);
}


