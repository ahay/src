/* Missing data interpolation in 2-D by Laplacian regularization. */
/*
  Copyright (C) 2014 University of Texas at Austin

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
#include "shapefill.h"

int main(int argc, char* argv[])
{
    int i, id, i3, n3, dim, dim1, n[SF_MAX_DIM], nd, rect[SF_MAX_DIM], niter;
    float *map, *msk;
    bool *known, verb;
    char key[6];
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("niter",&niter)) niter=200;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    
    dim = sf_filedims (in,n);
    if (!sf_getint("dim",&dim1)) dim1=dim;
    /* dimensionality */

    nd = 1;
    n3 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/
	if (i < dim1) {
	    nd *= n[i];
	} else {
	    n3 *= n[i];
	}
    }

    map = sf_floatalloc(nd);
    known = sf_boolalloc(nd);

    if (NULL != sf_getstring("mask")) {
	/* optional mask file with zeroes for missing data locations */
	mask = sf_input("mask");
	msk =  sf_floatalloc(nd);
    } else {
	mask = NULL;
	msk = map;
    }

    shapefill_init (dim, nd, n, rect, verb);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(map,nd,in);
	
	if (NULL != mask) sf_floatread(msk,nd,mask);

	for (id=0; id < nd; id++) {
	    known[id] = (bool) (msk[id] != 0.);
	}

	shapefill(niter,map,known);

	sf_floatwrite(map,nd,out);
    }

    exit(0);
}
