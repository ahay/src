/* Smooth derivative by plane-wave construction */
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

#include "smoothpwd.h"
#include "predict.h"


static int n;
static float *w, *p, *t;

static void predict_smooth_lop(bool adj, bool add, 
			       int nx, int ny, float* x, float* y);

void smoothpwd_init(int n1, int n2 /* data size */,
		    float eps      /* PWD regularization */,
		    int order      /* accuracy order */,
		    int rect1      /* vertical smoothing radius */, 
		    float **dip    /* dip field [n2][n1] */)
/*< initialize >*/
{
    n = n1*n2;

    sf_repeat_init(n1,n2,sf_causint_lop);
    predict_init (n1,n2,eps,order,1,false);
    predict_set(dip);
    sf_triangle2_init (rect1,1,n1,n2,1);

    w = sf_floatalloc(n);
    p = sf_floatalloc(n);
    t = sf_floatalloc(n);
}

void smoothpwd_close(void)
/*< free allocated storage >*/
{
    predict_close();
    free(w);
    free(p);
    free(t);
}

void smoothpwd(int niter     /* number of iterations */, 
	       int ncycle    /* number of cycles */, 
	       float* weight /* data weighting */, 
	       float* data   /* input data */, 
	       float* der    /* output derivative */,
	       bool verb     /* verbosity flag */,
	       float eps     /* regularization parameter */) 
/*< find the derivative >*/
{
    int i, iter;

    for (i=0; i < n; i++) {
	w[i] = 1.;
    }
    
    for (iter=0; iter < ncycle; iter++) {
	if (NULL != weight) {
	    sf_solver_prec (sf_repeat_lop,sf_cgstep,predict_smooth_lop,
			    n,n,n,der,data,niter,eps,
			    "wt",weight,"verb",verb,"mwt",w,"xp",p,"end");
	} else {
	    sf_solver_prec (sf_repeat_lop,sf_cgstep,predict_smooth_lop,
			    n,n,n,der,data,niter,eps,
			    "verb",verb,"mwt",w,"xp",p,"end");
	}
	sf_cgstep_close();
	
	for (i=0; i < n; i++) {
	    w[i] = p[i]; /* "Cauchy" weight + positivity */
	}	    
    }
}

static void predict_smooth_lop(bool adj, bool add, 
			       int nx, int ny, float* x, float* y)
{
    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	sf_triangle2_lop (true, false, nx, nx, t, y);
	predict_lop (true, add, nx, nx, x, t);
    } else {
	predict_lop (false, false, nx, nx, x, t);
	sf_triangle2_lop (false, add, nx, nx, t, y);
    }
}
