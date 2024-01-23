/* Random noise removal by nonstationary deconvolution */
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
#include "mmmult.h"

int main(int argc, char* argv[])
{
    int i, n1, n2, n12, n3, nf1, nf2, nf3, nf4, nf1234, i3, niter;
    float *data, *model, *filt, *weight, eps;
    sf_file inp, out, fil, wht;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    /* input data, output model */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(inp,2);

    /* input filter */
    if (!sf_histint(fil,"n1",&nf1)) sf_error("No n1= in filter");
    if (!sf_histint(fil,"n2",&nf2)) sf_error("No n2= in filter");
    if (!sf_histint(fil,"n3",&nf3)) sf_error("No n3= in filter");
    if (!sf_histint(fil,"n4",&nf4)) sf_error("No n4= in filter");
    
    if (nf3!=n1 || nf4!=n2) sf_error("need n1==nf3 && n2==nf4");
    nf1234 = nf1*nf2*nf3*nf4;
    
    filt = sf_floatalloc(nf1234);
    data = sf_floatalloc(n12);
    model = sf_floatalloc(n12);

    if (NULL != sf_getstring("weight")) {
	wht = sf_input("weight");
	weight = sf_floatalloc(n12);
	sf_weight_init(weight);
    } else {
	wht = NULL;
	weight = NULL;
    }

    mmmult_init (filt, nf1, nf2, nf3, nf4);	

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(filt,nf1234,fil);
	sf_floatread(data,n12,inp);

	if (NULL != wht) {
	    sf_floatread(weight,n12,wht);
	    for (i=0; i < n12; i++) {
		data[i] *= weight[i];
	    }
	    sf_solver_reg(sf_weight_lop,sf_cgstep,mmmult_lop,n12,n12,n12,
			  model,data,niter,1.0f,"verb",true,"end");
	} else {
	    sf_solver_reg(sf_copy_lop,sf_cgstep,mmmult_lop,n12,n12,n12,
			  model,data,niter,eps,"verb",true,"end");
	}
	sf_cgstep_close();

	sf_floatwrite(model,n12,out);
    }
    
    exit(0);
}

