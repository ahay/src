/* Missing data interpolation using one or two prediction-error filters. */
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

#include "createhelix.h"
#include "bound.h"
#include "misinput.h"
#include "pef.h"
#include "compress.h"
#include "printfilter.h"
#include "mask2i.h"

int main(int argc, char* argv[])
{
    int dim, i, niter, n12, n[SF_MAX_DIM], a[SF_MAX_DIM], b[SF_MAX_DIM];
    int center[SF_MAX_DIM], gap[SF_MAX_DIM], *known;
    float *dd, *hh;
    sf_filter aa, bb;
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_input("mask");
    if (SF_INT != sf_gettype(mask)) sf_error("Need int type in mask");

    if (!sf_getint("niter",&niter)) niter=80;
    /* number of iterations */

    dim = sf_filedims (in,n);
    n12 = 1;
    for( i=0; i < dim; i++) {
	n12 *= n[i];
	gap[i] = 0;
    }
    if (!sf_getints("center",center,dim)) {
	/* filter center */
	for( i=0; i < dim; i++) {
	    center[i]=0;
	}
    }
  
    dd = sf_floatalloc(n12);
    hh = sf_floatalloc(n12);
    known = sf_intalloc(n12);

    sf_floatread(dd,n12,in);
    sf_intread(known,n12,mask);

    if (!sf_getints("a",a,dim)) sf_error("Need a=");
    /* first filter dimensions */

    aa = createhelix(dim,n,center,gap,a);

    bound(dim, n, n, a, aa);
    find_mask(n12, known, aa);
    find_pef (n12, dd, aa, 2*aa->nh);
    aa = compress(aa, FLT_EPSILON);
    print(dim, n, center, a, aa);

    if (!sf_getints("b",b,dim)) {
    /* second filter dimensions */
	bb = NULL;
    } else {
	bb = createhelix(dim,n,center,gap,b);

	bound(dim, n, n, b, bb);
	find_mask(n12, known, bb);
	find_pef (n12, dd, bb, 2*bb->nh);
	bb = compress(bb, FLT_EPSILON);
	print(dim, n, center, b, bb);
    }

    maski (niter, n12, dd, hh, known, aa, bb);
    sf_floatwrite(hh,n12,out);

    exit(0);
}

