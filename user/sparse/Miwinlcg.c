/* Image-domain waveform tomography (Non-linear CG). */
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

#include "iwinlcg.h"

int main(int argc, char* argv[])
{
    bool verb, load;
    int n1, n2, npml, nh, ns, nw;
    int prect[3], pliter;
    int dorder, grect[2], gliter;
    float plower, pupper;
    float vpml, d1, d2, **vel, dw, ow;
    float ***mask, ***wght, **prec;
    char *datapath;
    sf_file in, out, source, data;
    sf_file imask, weight, precon;
    int uts, mts, i, j;
    char *order;
    float fx0, fx1, fx2, *x0, *x1, *g0, *g1, *s, beta;
    float geps, gscale, delta, lower, upper, alpha, ltol, armijo;
    int iter, miter, liter;

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

    if (!sf_getfloat("vpml",&vpml)) vpml=4.;
    /* PML velocity */

    if (!sf_getint("npml",&npml)) npml=20;
    /* PML width */

    if (NULL == (order = sf_getstring("order"))) order="j";
    /* discretization scheme (default optimal 9-point) */

    if (!sf_getint("prect1",&prect[0])) prect[0]=5;
    /* slope smoothing radius on axis 1 */
    if (!sf_getint("prect2",&prect[1])) prect[1]=1;
    /* slope smoothing radius on axis 2 */
    if (!sf_getint("prect3",&prect[2])) prect[2]=5;
    /* slope smoothing radius on axis 3 */

    if (!sf_getint("pliter",&pliter)) pliter=20;
    /* slope estimation # of linear iterations */

    if (!sf_getfloat("plower",&plower)) plower=0.1;
    /* slope thresholding lower limit */
    if (!sf_getfloat("pupper",&pupper)) pupper=4.;
    /* slope thresholding upper limit */

    if (!sf_getint("dorder",&dorder)) dorder=6;
    /* image derivative accuracy order */
    
    if (!sf_getint("grect1",&grect[0])) grect[0]=5;
    /* gradient smoothing radius on axis 1 */
    if (!sf_getint("grect2",&grect[1])) grect[1]=5;
    /* gradient smoothing radius on axis 2 */

    if (!sf_getint("gliter",&gliter)) gliter=1;
    /* # of Gauss-Newton iterations */

    if (!sf_getfloat("geps",&geps)) geps=0.;
    /* regularization parameter for Gauss-Newton */

    if (!sf_getfloat("gscale",&gscale)) gscale=0.5;
    /* gradient re-scale */

    if (!sf_getfloat("delta",&delta)) delta=0.05;
    /* Nonlinear-CG termination */

    if (!sf_getfloat("ltol",&ltol)) ltol=1.e-4;
    /* Armijo condition parameter */

    if (!sf_getint("miter",&miter)) miter=10;
    /* Nonlinear-CG maximum # of iterations */

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
    if (NULL == sf_getstring("weight")) {
	weight = NULL;
	wght = NULL;
    } else {
	weight = sf_input("weight");
	wght = sf_floatalloc3(n1,n2,2*nh+1);
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
    x0 = sf_floatalloc(n1*n2);
    x1 = sf_floatalloc(n1*n2);
    g0 = sf_floatalloc(n1*n2);
    g1 = sf_floatalloc(n1*n2);
    s  = sf_floatalloc(n1*n2);

    for (j=0; j < n2; j++) {
	for (i=0; i < n1; i++) {
	    x0[j*n1+i] = vel[j][i];
	}
    }

    /* initialize operators */
    iwinlcg_init(false,order, npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 source,data, load,datapath, uts,
		 prect[0],prect[1],prect[2],
		 pliter,plower,pupper,
		 dorder,
		 grect[0],grect[1],
		 gliter,geps,gscale,
		 lower,upper);
    
    /* non-linear CG optimization */
    fx0 = iwinlcg_eval(x0, mask,wght);
    iwinlcg_grad(x0, wght,prec,g0);

    for (j=0; j < n2; j++) {
	for (i=0; i < n1; i++) {
	    s[j*n1+i] = -g0[j*n1+i];
	}
    }

    for (iter=0; iter < miter; iter++) {
	if (verb) sf_warning("Iteration %d...",iter+1);

	/* line-search */
	alpha = 1.; fx2 = SF_HUGE;
	for (liter=0; liter < 10; liter++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    x1[j*n1+i] = x0[j*n1+i]+alpha*s[j*n1+i];
		}
	    }

	    fx1 = iwinlcg_eval(x1, mask,wght);
	    iwinlcg_grad(x1, wght,prec,g1);

	    armijo = alpha*cblas_sdot(n1*n2,s,1,g0,1);

	    if (verb) sf_warning("Line search %d: fx0=%g, fx1=%g, armijo=%g."
				 ,liter+1,fx0,fx1,armijo);

	    if (fx1 <= fx0+ltol*armijo)
		break;
	    
	    if (fx1 < fx2) {
		alpha *= 0.5;
	    } else {
		liter = 10;
		break;
	    }

	    fx2 = fx1;
	}
	if (liter == 10) {
	    sf_warning("Iteration terminated due to line search failure.");
	    
	    sf_floatwrite(x0,n1*n2,out);
	    iwinlcg_free();

	    exit(0);
	}

	/* replace model */
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		x0[j*n1+i] = x1[j*n1+i];
	    }
	}

	/* test convergence */
	if (cblas_snrm2(n1*n2,g1,1) <= 
	    delta*SF_MAX(1.,cblas_snrm2(n1*n2,x1,1))) {
	    if (verb) sf_warning("Iteration terminated due to convergence.");
	    break;
	}

	beta = SF_MAX(0.,(cblas_sdot(n1*n2,g1,1,g1,1)
			  -cblas_sdot(n1*n2,g1,1,g0,1))
		      /cblas_sdot(n1*n2,g0,1,g0,1));

	if (verb) sf_warning("Done with update beta=%g.",beta);

	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		s[j*n1+i] = -g1[j*n1+i]+beta*s[j*n1+i];
		g0[j*n1+i] = g1[j*n1+i];
	    }
	}
	fx0 = fx1;
    }

    sf_floatwrite(x0,n1*n2,out);
    iwinlcg_free();

    exit(0);
}
