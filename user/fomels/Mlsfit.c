/* Simple least-squares regression. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "gaussel.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], dim, dim1, n1, n2, i, i2, i1, ic, id, nc;
    float *dat, **func, **wfunc, **mat, *rhs, *sol, *weight;
    sf_file inp, fit, coef, out, wht;

    sf_init(argc,argv);
    inp = sf_input("in");
    fit = sf_input("fit");
    out = sf_output("out");

    dim = sf_filedims(inp,n);
    if (!sf_getint("dim",&dim1)) dim1=dim;
    /* number of dimensions */
 
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i < dim1) 
	    n1 *= n[i];
	else 
	    n2 *= n[i];
    }

    nc = sf_filesize(fit)/n1;

    if (NULL != sf_getstring("coef")) {
	coef = sf_output("coef");
	sf_putint(coef,"n1",nc);
	for (i=1; i < dim1; i++) {
	    sf_unshiftdim(coef,coef,2);
	}
    } else {
	coef = NULL;
    }

    dat = sf_floatalloc(n1);
    func = sf_floatalloc2(n1,nc);

    sf_floatread(func[0],n1*nc,fit);
    sf_fileclose(fit);

    mat = sf_floatalloc2(nc,nc);
    
    if (NULL != sf_getstring("weight")) {
      wht = sf_input("weight");
      weight = sf_floatalloc(n1);
      wfunc = sf_floatalloc2(n1,nc);
    } else {
      wht = NULL;
      weight = NULL;
      wfunc = func;

	/* compute A'A matrix */
	for (ic=0; ic < nc; ic++) {
	  for (id=0; id <= ic; id++) {
	    mat[ic][id] = 0.;
	    for (i1=0; i1 < n1; i1++) {
	      mat[ic][id] += func[ic][i1]*func[id][i1];
	    }
	    mat[id][ic] = mat[ic][id];
	  }
	}
    }

    gaussel_init(nc);
    sol = sf_floatalloc(nc);
    rhs = sf_floatalloc(nc);
    mat = sf_floatalloc2(nc,nc);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(dat,n1,inp);
	
	if (NULL != weight) {
	  sf_floatread(weight,n1,wht);

	  for (ic=0; ic < nc; ic++) {
	    for (i1=0; i1 < n1; i1++) {
	      wfunc[ic][i1] = func[ic][i1]*weight[i1];
	    }
	  }

	  /* compute A'A matrix */
	  for (ic=0; ic < nc; ic++) {
	    for (id=0; id <= ic; id++) {
	      mat[ic][id] = 0.;
	      for (i1=0; i1 < n1; i1++) {
		mat[ic][id] += wfunc[ic][i1]*wfunc[id][i1];
	      }
	      mat[id][ic] = mat[ic][id];
	    }
	  }
	  
	  for (i1=0; i1 < n1; i1++) {
	    dat[i1] *= weight[i1];
	  }
	}	

	/* compute A'd */
	for (ic=0; ic < nc; ic++) {
	    rhs[ic] = 0.;
	    for (i1=0; i1 < n1; i1++) {
		rhs[ic] += wfunc[ic][i1]*dat[i1];
	    }
	}

	/* inversion */
	gaussel_solve(mat,rhs,sol);

	if (NULL != coef) sf_floatwrite(sol,nc,coef);

	/* compute Ac */
	for (i1=0; i1 < n1; i1++) {
	    dat[i1] = 0.;
	    for (ic=0; ic < nc; ic++) {
	      dat[i1] += func[ic][i1]*sol[ic];
	    }
	}

	sf_floatwrite(dat,n1,out);
    }

    exit(0);
}
