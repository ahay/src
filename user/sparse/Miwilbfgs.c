/* Image-domain waveform tomography (L-BFGS). */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lbfgs.h"
#include "iwilbfgs.h"

static float ***mask, ***wght, **prec;

static lbfgsfloatval_t evaluate(void *instance,
				const lbfgsfloatval_t *x,
				lbfgsfloatval_t *g,
				const int n,
				const lbfgsfloatval_t step)
/* evaluate objective function and gradient */
{
    lbfgsfloatval_t fx;

    fx = iwilbfgs_eval(x,mask,wght);

    iwilbfgs_grad(x,wght,prec,g);

    return fx;
}

static int progress(void *instance,
		    const lbfgsfloatval_t *x,
		    const lbfgsfloatval_t *g,
		    const lbfgsfloatval_t fx,
		    const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm,
		    const lbfgsfloatval_t step,
		    int n,
		    int k,
		    int ls)
/* report optimization progress */
{
    sf_warning("L-BFGS iteration %d: fx=%g, xnorm=%g, gnorm=%g, step=%g after %d."
	       ,k,fx,xnorm,gnorm,step,ls);

    return 0;
}

int main(int argc, char* argv[])
{
    bool verb, load;
    int n1, n2, npml, nh, ns, nw;
    int grect[2], mline;
    float gscale, epsilon;
    float d1, d2, **vel, dw, ow;
    char *datapath;
    sf_file in, out, source, data;
    sf_file imask, weight, precon;
    int uts, mts, ret, i, j;
    char *order;
    lbfgsfloatval_t fx, *x;
    lbfgs_parameter_t param;
    float lower, upper;
    int nhess, miter;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("nh",&nh)) nh=0;
    /* horizontal space-lag */

    if (!sf_getbool("load",&load)) load=false;
    /* load LU */

    if (!sf_getint("uts",&uts)) uts=0;
    /* number of OMP threads */

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    uts = (uts < 1)? mts: uts;

    if (!sf_getint("npml",&npml)) npml=10;
    /* PML width */

    if (NULL == (order = sf_getstring("order"))) order="j";
    /* discretization scheme (default optimal 9-point) */    

    if (!sf_getint("grect1",&grect[0])) grect[0]=5;
    /* gradient smoothing radius on axis 1 */
    if (!sf_getint("grect2",&grect[1])) grect[1]=5;
    /* gradient smoothing radius on axis 2 */

    if (!sf_getfloat("gscale",&gscale)) gscale=0.5;
    /* gradient re-scale */

    if (!sf_getint("nhess",&nhess)) nhess=6;
    /* L-BFGS # of Hessian corrections */

    if (!sf_getint("miter",&miter)) miter=10;
    /* L-BFGS maximum # of iterations */

    if (!sf_getint("mline",&mline)) mline=5;
    /* L-BFGS maximum # of line search */

    if (!sf_getfloat("epsilon",&epsilon)) epsilon=1.e-2;
    /* L-BFGS termination epsilon */

    if (!sf_getfloat("lower",&lower)) lower=1.5;
    /* lower bound of feasible set */

    if (!sf_getfloat("upper",&upper)) upper=7.5;
    /* upper bound of feasible set */

    /* read initial model */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

    vel = sf_floatalloc2(n1,n2);
    sf_floatread(vel[0],n1*n2,in);

    if (load)
	datapath = sf_histstring(in,"in");
    else
	datapath = NULL;

    /* read source */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histint(source,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(source,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(source,"o4",&ow)) sf_error("No ow=.");

    /* read data */
    if (NULL == sf_getstring("data"))
	sf_error("Need data=");
    data = sf_input("data");

    /* read image mask */
    if (NULL == sf_getstring("imask")) {
	imask = NULL;
	mask = NULL;
    } else {
	imask = sf_input("imask");
	mask = sf_floatalloc3(n1,n2,2*nh+1);
	sf_floatread(mask[0][0],n1*n2*(2*nh+1),imask);
	sf_fileclose(imask);
    }

    /* read image weight */
    wght = sf_floatalloc3(n1,n2,2*nh+1);    

    if (NULL == sf_getstring("weight")) {
	weight = NULL;
	for (i=0; i < n1*n2*(2*nh+1); i++) {
	    wght[0][0][i] = 1.;
	}
    } else {
	weight = sf_input("weight");
	sf_floatread(wght[0][0],n1*n2*(2*nh+1),weight);
	sf_fileclose(weight);
    }

    /* read model preconditioner */
    if (NULL == sf_getstring("precon")) {
	precon = NULL;
	prec = NULL;
    } else {
	precon = sf_input("precon");
	prec = sf_floatalloc2(n1,n2);
	sf_floatread(prec[0],n1*n2,precon);
	sf_fileclose(precon);
    }

    /* allocate temporary memory */
    x = lbfgs_malloc(n1*n2);

    for (j=0; j < n2; j++) {
	for (i=0; i < n1; i++) {
	    x[j*n1+i] = (lbfgsfloatval_t) vel[j][i];
	}
    }

    /* initialize operators */
    iwilbfgs_init(verb,order, npml,
		  n1,n2, d1,d2,
		  nh,ns, ow,dw,nw,
		  source,data, load,datapath, uts,
		  grect[0],grect[1],
		  gscale,
		  lower,upper);
    
    /* initialize L-BFGS */
    lbfgs_parameter_init(&param);

    param.m = nhess;
    param.max_iterations = miter;
    param.max_linesearch = mline;
    param.epsilon = (lbfgsfloatval_t) epsilon;

    /* L-BFGS optimization */
    if (verb) {
	ret = lbfgs(n1*n2,x,&fx, evaluate,progress, NULL,&param);
    } else {
	ret = lbfgs(n1*n2,x,&fx, evaluate,NULL,     NULL,&param);
    }

    if (verb) sf_warning("L-BFGS optimization terminated with status code %d.",ret);

    /* write output */
    for (j=0; j < n2; j++) {
	for (i=0; i < n1; i++) {
	    vel[j][i] = (float) x[j*n1+i];
	}
    }

    sf_floatwrite(vel[0],n1*n2,out);

    /* free */
    iwilbfgs_free();

    exit(0);
}
