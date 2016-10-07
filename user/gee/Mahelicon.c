/* Apply multidimensional nonstationary filter on a helix. */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

#include "regrid.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], m[SF_MAX_DIM];
    int ndim, dim, n123, n123s, i, ia, ns, i1, na, i4, n4;
    float *f, *dd;
    char *lagfile;
    sf_filter aa;
    sf_file in, filt, lag, out;
 
    sf_init(argc,argv);

    in = sf_input("in");
    filt = sf_output("filt");
    out = sf_output("out");

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");

    aa = sf_allocatehelix (na);

    if (NULL != (lagfile = sf_getstring("lag")) 
	/*( lag file with filter lags )*/
	|| 
	NULL != (lagfile = sf_histstring(filt,"lag"))) {
	lag = sf_input(lagfile);

	sf_intread(aa->lag,na,lag);
    } else {
	lag = NULL;
	for( ia=0; ia < na; ia++) {
	    aa->lag[ia] = ia+1;
	}
    }


    ndim = sf_filedims(in,n);

    if (!sf_getint("dim",&dim)) dim=ndim; /* number of dimensions */
    
    if (!sf_getints ("n",m,dim) && (NULL == lag ||
				    !sf_histints (lag,"n",m,dim))) {
	for (i=0; i < dim; i++) {
	    m[i] = n[i];
	}
    }

    regrid (dim, m, n, aa);

    n4 = sf_leftsize(in,dim);
    
    n123 = 1;
    for (i=0; i < dim; i++) {
	n123 *= n[i];
    }

    dd = sf_floatalloc(n123);

    n123s = n123*na;    
    f = sf_floatalloc(n123s);

    for (i4=0; i4 < n4; i4++) {

	sf_floatread(dd,n123,in);
	sf_floatread(f,n123s,in);
	
	/* apply shifts: dd -> d */
	for (i=ia=0; ia < na; ia++) {
	    ns = aa->lag[ia];
	    for (i1=0; i1 < n123; i1++,i++) {
		dd[i1] -= f[i]*dd[i1-ns];
	    }
	}

	sf_floatwrite(dd,n123,out);
    }
    
    exit(0);
}
