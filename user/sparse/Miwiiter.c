/* Image-domain waveform tomography. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "iwioper.h"

int main(int argc, char* argv[])
{
    bool verb, load, mass, shape;
    int n1, n2, cgiter, npml; 
    int nh, ns, nw, n[2], rect[2];
    float d1, d2, **vel, dw, ow, tol;
    float *di, *dm, ***wght, reg, **prec, **xmov, **rmov, *p=NULL;
    char *datapath;
    sf_file in, out, model, us, ur;
    sf_file weight, precon, miter, riter;
    int uts, mts;
    char *order;

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

    if (!sf_getbool("mass",&mass)) mass=false;
    /* if y, use discretization-based mass term */

    if (!sf_getint("cgiter",&cgiter)) cgiter=10;
    /* number of conjugate-gradient iterations */

    if (!sf_getbool("shape",&shape)) shape=false;
    /* regularization (default Tikhnov) */

    if (!sf_getfloat("reg",&reg)) reg=0.;
    /* regularization parameter */

    /* read model */
    if (NULL == sf_getstring("model"))
	sf_error("Need model=");
    model = sf_input("model");

    if (!sf_histint(model,"n1",&n1)) sf_error("No n1= in model.");
    if (!sf_histint(model,"n2",&n2)) sf_error("No n2= in model.");

    if (!sf_histfloat(model,"d1",&d1)) sf_error("No d1= in model.");
    if (!sf_histfloat(model,"d2",&d2)) sf_error("No d2= in model.");

    vel = sf_floatalloc2(n1,n2);
    sf_floatread(vel[0],n1*n2,model);

    if (load)
	datapath = sf_histstring(model,"in");	
    else
	datapath = NULL;

    /* allocate memory */
    dm = sf_floatalloc(n1*n2);
    di = sf_floatalloc(n1*n2*(2*nh+1));

    /* read input */
    sf_floatread(di,n1*n2*(2*nh+1),in);

    /* read source wavefield */
    if (NULL == sf_getstring("us"))
	sf_error("Need source wavefield us=");
    us = sf_input("us");

    if (!sf_histint(us,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histint(us,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(us,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(us,"o4",&ow)) sf_error("No ow=.");    
    
    /* read receiver wavefield */
    if (NULL == sf_getstring("ur"))
	sf_error("Need receiver wavefield ur=");
    ur = sf_input("ur");

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

    /* write output header */
    sf_putint(out,"n3",1);

    /* output model and residual iteration */
    if (NULL == sf_getstring("miter") && 
	NULL == sf_getstring("riter")) {
	miter = NULL;
	riter = NULL;
	xmov = NULL;
	rmov = NULL;
    } else {
	miter = sf_output("miter");
	sf_putint(miter,"n3",cgiter);
	riter = sf_output("riter");
	sf_putint(riter,"n4",cgiter);
	xmov = sf_floatalloc2(n1*n2,cgiter);
	rmov = sf_floatalloc2(n1*n2*(2*nh+1),cgiter);
    }    

    /* initialize operator */
    iwi_init(npml, n1,n2,d1,d2, nh,ns,ow,dw,nw,
	     us,ur, load,datapath, uts, order,mass);

    /* initialize regularization */
    if (shape) {
	n[0] = n1; n[1] = n2;

	if (!sf_getfloat("tol",&tol)) tol=1.e-6;
	/* tolerance for shaping regularization */

	if (!sf_getint("rect1",&rect[0])) rect[0]=1;
	/* smoothing radius on axis 1 */
	if (!sf_getint("rect2",&rect[1])) rect[1]=1;
	/* smoothing radius on axis 2 */
	
	/* triangle smoothing operator */
	sf_trianglen_init(2,rect,n);
	sf_repeat_init(n1*n2,1,sf_trianglen_lop);

	sf_conjgrad_init(n1*n2,n1*n2,n1*n2*(2*nh+1),n1*n2*(2*nh+1),reg,tol,verb,false);
	p = sf_floatalloc(n1*n2);
    } else {
	/* 2D gradient operator */
	sf_igrad2_init(n1,n2);
    }

    /* set weight and preconditioner */
    iwi_set(vel,wght,prec);

    /* solve update */
    if (shape) {
	sf_conjgrad(NULL,iwi_oper,sf_repeat_lop,p,dm,di,cgiter);
    } else {
	if (NULL == miter && NULL == riter) {
	    sf_solver_reg(iwi_oper,sf_cgstep,sf_igrad2_lop,
			  2*n1*n2,n1*n2,n1*n2*(2*nh+1),dm,di,
			  cgiter,reg,"verb",verb,"end");
	} else {
	    sf_solver_reg(iwi_oper,sf_cgstep,sf_igrad2_lop,
			  2*n1*n2,n1*n2,n1*n2*(2*nh+1),dm,di,
			  cgiter,reg,"xmov",xmov,"rmov",rmov,"verb",verb,"end");
	}

	sf_cgstep_close();
    }
    
    /* write output */
    sf_floatwrite(dm,n1*n2,out);

    if (NULL != miter || NULL != riter) {
	sf_floatwrite(xmov[0],n1*n2*cgiter,miter);
	sf_floatwrite(rmov[0],n1*n2*(2*nh+1)*cgiter,riter);
    }

    exit(0);
}
