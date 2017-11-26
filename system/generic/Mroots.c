/* Find roots of a complex polynomial. */
/*
  Copyright (C) 2007 The University of Texas at Austin

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
#include "ceig.h"

static float (*func)(sf_complex);

static int compare(const void * a, const void * b)
{
    float aa, bb;

    aa = func(* (sf_complex*) a);
    bb = func(* (sf_complex*) b);

    return (aa<bb)? -1: (aa>bb)? 1:0;
}

int main(int argc, char* argv[])
{
    bool verb;
    char *sort;
    int j, k, n, m, i2, n2, niter, *map;
    sf_complex **a, *e, *old;
    float tol, dist, dk;
    sf_file poly, root;

    sf_init(argc,argv);
    poly = sf_input("in");
    root = sf_output("out");

    if (SF_COMPLEX != sf_gettype(poly)) sf_error("Need complex input");
    if (!sf_histint(poly,"n1",&n)) sf_error("No n1= in input");
    n2 = sf_leftsize(poly,1);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    if (!sf_getfloat("tol",&tol)) tol=1.0e-6;
    /* tolerance for convergence */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (NULL == (sort = sf_getstring("sort"))) sort="real";
    /* attribute for sorting (phase,amplitude,real,imaginary) */

    switch (sort[0]) {
	case 'p':
	    func = cargf;
	    break;
	case 'a':
	    func = cabsf;
	    break;
	case 'i':
	    func = cimagf;
	    break;
	case 'r':
	default:
	    func = crealf;
	break;
    }

    sf_putint(root,"n1",n-1);

    a = sf_complexalloc2(n,n);
    e = sf_complexalloc(n);
    old = sf_complexalloc(n-1);
    map = sf_intalloc(n-1);

    ceig_init(verb,n);

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(e,n,poly);

	for (m = n; m > 0; m--) {
	    if (cabsf(e[m-1]) > FLT_EPSILON) break;
	}
	m--;

	for (j=0; j < m; j++) {
	    for (k=0; k < m; k++) {
		a[j][k] = sf_cmplx(0.,0.);
	    }
	}
	for (j=0; j < m-1; j++) {
	    a[j][j+1]=sf_cmplx(1.,0.);
	}
	for (j=0; j < m; j++) {
#ifdef SF_HAS_COMPLEX_H
	    a[m-1][j]=-e[j]/e[m];
#else
	    a[m-1][j]=sf_cneg(sf_cdiv(e[j],e[m]));
#endif
	}

	ceig(niter,tol,m,a,e);
	
	if (0==i2) {
	    /* sort the first set of roots */
	    qsort(e,n-1,sizeof(sf_complex),compare);
	    for (j=0; j < n-1; j++) {
		old[j]=e[j];
	    }
	} else {
	    /* find nearest to previous */
	    for (j=0; j < n-1; j++) {
		map[j] = -1;
	    }
	    /* loop through old roots */
	    for (j=0; j < n-1; j++) {
		dist = SF_HUGE;
		/* find nearest not taken */
		for (k=0; k < n-1; k++) {
		    if (map[k] >= 0) continue;
		    /* Euclidean distance */
#ifdef SF_HAS_COMPLEX_H
		    dk = cabsf(old[j]-e[k]);
#else
		    dk = cabsf(sf_cadd(old[j],sf_crmul(e[k],-1.0)));
#endif
		    if (dk < dist) {
			m = k;
			dist = dk;
		    }
		}
		map[m] = j;
		old[j] = e[m];
	    }
	}

	sf_complexwrite(old,n-1, root);
    }

    exit(0);
}
