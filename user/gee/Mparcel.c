/* Patching test.
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

#include <rsf.h>

#include "parcel.h"

int main (int argc, char* argv[]) {
    int dim, i, j, n[SF_MAX_DIM], w[SF_MAX_DIM], k[SF_MAX_DIM], np, nw, n12;
    float *data, *wind;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    dim = sf_filedims(in,n);
    if (!sf_getints("w",w,dim)) sf_error("Need w=");
    if (!sf_getints("k",k,dim)) sf_error("Need k=");

    n12 = 1; np = 1; nw = 1;
    for(j=0; j < dim; j++) {
	n12 *= n[j];
	np *= k[j];
	nw *= w[j];
    }
    nw *= np;

    data = sf_floatalloc (n12);
    wind = sf_floatalloc (nw);

    sf_floatread(data,n12,in);

    parcel_init (dim, k, n, w);
    parcel_lop (false, false, n12, nw, data, wind);
    parcel_lop ( true, false, n12, nw, data, wind);

    for (i=0; i < n12; i++) {
	sf_line2cart(dim, n, i, k);
	if (k[0] <= 1 || k[0] >= n[0] - 2) data[i] = 0.;
    }

    sf_floatwrite (data,n12,out);

    exit (0);
}

/* 	$Id: Mparcel.c 4796 2009-09-29 19:39:07Z ivlad $	 */
