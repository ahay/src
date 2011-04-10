/* Moveout flattening. */
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
#include "predict.h"

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, i3, i0, order;
    float eps, **u, **p, **v;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("i0",&i0)) i0=0;
    /* reference trace */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1, n2, eps*eps, order, 1);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    v = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	if (verb) fprintf(stderr,"cmp %d of %d\n",i3+1,n3);

	sf_floatread(u[0],n1*n2,in);
	sf_floatread(p[0],n1*n2,dip);

	predict_flat(i0, u, v, p);

	sf_floatwrite(v[0],n1*n2,out);
    }

    exit (0);
}

/* 	$Id$	 */
