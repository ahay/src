/* Random noise removal by deconvolution on a helix */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "regrid.h"

int main(int argc, char* argv[])
{
    int dim, na, nx, i, j, niter, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float *data, *model, *weight, eps;
    char *lagfile;
    sf_filter pef;
    sf_file inp, out, fil, lag, wht;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    dim = sf_filedims(inp,n);
    
    if (!sf_histint(fil,"n1",&na)) sf_error("No n1= in sfilt");

    pef = sf_allocatehelix(na);

    if (NULL == (lagfile = sf_histstring(fil,"lag")) &&
	NULL == (lagfile = sf_getstring("lag"))) 
	sf_error("Need lag=");
    lag = sf_input(lagfile);
    if (!sf_histints(lag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(pef->lag,na,lag);
    regrid(dim,m,n,pef);
    sf_fileclose(lag);

    sf_helicon_init (pef);
    
    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    /* input data, output model */
    nx=1;
    for(j=0; j < dim; j++) {
	nx *= n[j];
    }

    data = sf_floatalloc(nx);
    model = sf_floatalloc(nx);

    if (NULL != sf_getstring("weight")) {
	wht = sf_input("weight");
	weight = sf_floatalloc(nx);
	sf_weight_init(weight);
    } else {
	wht = NULL;
	weight = NULL;
    }


    sf_floatread(pef->flt,na,fil);
    sf_floatread(data,nx,inp);

    if (NULL != wht) {
	sf_floatread(weight,nx,wht);
	for (i=0; i < nx; i++) {
	    data[i] *= weight[i];
	}
	sf_solver_reg(sf_weight_lop,sf_cgstep,sf_helicon_lop,nx,nx,nx,
		      model,data,niter,1.0f,"verb",true,"end");
    } else {
	sf_solver_reg(sf_copy_lop,sf_cgstep,sf_helicon_lop,nx,nx,nx,
		      model,data,niter,eps,"verb",true,"end");
    }

    sf_floatwrite(model,nx,out);
    
    exit(0);
}

