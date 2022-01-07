/* Diffusion by regularized inversion. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <assert.h>

#include <rsf.h>

static int n1, n2, n12;

static void gradient (bool adj, bool add, 
		      int nx, int ng, float* x, float* g)
/*< gradient operator >*/
{
    int i1,i2, i;
    assert(n12 == n1*n2 && nx == n12 && ng == 2*n12);

    sf_adjnull (adj,add,nx,ng,x,g);
    
    for (i2=1; i2 < n2-1; i2++) {  
	for (i1=1; i1 < n1-1; i1++) {
	    i = i1+i2*n1; /* map 2-D to 1-D */

	    if (adj) {
	      x[i+1] += 0.5*g[i];
	      x[i-1] -= 0.5*g[i];
	      x[i+n1] += 0.5*g[i+n12];
	      x[i-n1] -= 0.5*g[i+n12];
	    } else {
		g[i]     += 0.5*(x[i+1 ] - x[i-1]); 
		g[i+n12] += 0.5*(x[i+n1] - x[i-n1]);
	    }
	}
    }
}

int main(int argc, char* argv[])
{
    int niter, i, repeat;
    float *data, eps;
    double dot1[2], dot2[2];
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* Get dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2=");
    n12 = n1*n2;

    /* Allocate and read data */
    data = sf_floatalloc(n12);
    sf_floatread(data,n12,in);

    /* First let us do the dot-product test */
    sf_dot_test(gradient,n12,2*n12,dot1,dot2);
    sf_warning("Dot product test");
    sf_warning("****************");
    sf_warning("Check if %g == %g",dot1[0],dot1[1]);
    sf_warning("Check if %g == %g",dot2[0],dot2[1]);
    sf_warning("****************");
    
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of conjugate-gradient iterations */
    if (!sf_getint("repeat",&repeat)) repeat=1;
    /* number of smoothing iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */
    //return;
   for (i=0; i < repeat; i++) {
	/* Call a generic function for solving
	   .   F[m] =~ d
	   eps*D[m] =~ 0 
	*/
	sf_solver_reg(sf_copy_lop /* F */,
		      sf_cgstep   /* CG step */,
		      gradient    /* D */,
		      2*n12       /* size of D[m] */,
		      n12         /* size of m */,
		      n12         /* size of d */,
		      data        /* m (output) */,
		      data        /* d (input)  */,
		      niter       /* # iterations */,
		      eps         /* eps */,
		      "verb",true,"end");
	sf_cgstep_close(); /* free internal storage */
    }

    sf_floatwrite(data,n12,out);

    exit(0);
}
	
