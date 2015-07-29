/* Remove bursty noise by IRLS. */
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

#include "deburst.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2, niter;
    char *norm;
    float *data, *model, eps;
    sf_file in, out;
    sf_weight weight=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    model = sf_floatalloc(n1);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */
    if (NULL == (norm = sf_getstring("norm"))) {
		/* norm to use in IRLS (cauchy,l1) */
		weight=sf_cauchy;
	} else {
		sf_warning("got %s",norm);

		switch(norm[0]) {
			case 'c': case 'C':
				weight=sf_cauchy;
				break;
			case 'l': case 'L':
				weight=sf_l1;
				break;
			default:
				sf_error("unknown norm %s",norm);
				break;
		}
    }

    sf_irls_init(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (data,n1,in);

	deburst (n1, niter, weight, eps, data, model);

	sf_floatwrite (model,n1,out);
    }

    exit(0);
}

/* 	$Id: Mdeburst.c 7267 2011-06-13 18:38:18Z saragiotis $	 */
