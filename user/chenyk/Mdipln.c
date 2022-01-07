/* large dip calculation via non-stationary regularization  */
/*
  Copyright (C) 2012 The University of Texas at Austin 
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>
#include "dipln.h"

int main(int argc, char*argv[])
{
	sf_file in, out;
	int m, n, n1, n2, n3, niter, liter;
	int i3;
	bool verb, slope;
	float **wav, **dip, radius, eta, dip0;
	char *interp;

	/*added for non-stationary regularization*/
    int   *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */	
	int i, j, b, n123, dim;
	float eps;
    char key[8];
    int box[SF_MAX_DIM], nn[SF_MAX_DIM];
	sf_file rect[SF_MAX_DIM], shift[SF_MAX_DIM]; 
	/*added for non-stationary regularization*/
	
	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in, 2);

	if(!sf_getint("m", &m)) m=1;
	/* b[-m, ... ,n] */
	if(!sf_getint("n", &n)) n=1;
	/* b[-m, ... ,n] */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	/* nf = m+n+1; */
	
	if (!sf_getint("niter",&niter)) niter=5;
	/* number of iterations */
	if (!sf_getint("liter",&liter)) liter=20;
	/* number of linear iterations */
	if (!sf_getfloat("radius", &radius)) radius = 1.0;
	/* interpolating radius for opwd */
	if (!sf_getfloat("eta", &eta)) eta = 0.5;
	/* steps for iteration */
	if (!sf_getfloat("dip0", &dip0)) dip0 = 0.0;
	/* starting dip */
	if (!sf_getbool("verb", &verb)) verb = false;
	/* verbosity flag */
	if (!sf_getbool("slope", &slope)) slope = false;
	/* slope (y) or dip (n) estimation */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */


	wav = sf_floatalloc2(n1, n2);
	dip = sf_floatalloc2(n1, n2);

	/*added for non-stationary regularization*/
	dim = sf_filedims(in,nn);
    if (dim < 2) nn[1]=1;
    if (dim < 3) {nn[2]=1; box[2]=1;}
    n123 = nn[0]*nn[1]*nn[2];
    
	sf_warning("dim=%d, n123=%d",dim,n123);
	
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    rect[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rect[i] = NULL;
	    shift[i] = NULL;
	}
    }

	/*reading the non-stationary smoothing radii*/
    for (i=0; i < dim; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_floatalloc (n123);
	    sft[i] = sf_intalloc (n123);
		/* non-stationary dip smoothness on 1st/2nd/3rd axis */

	    sf_floatread(rct[i],n123,rect[i]);
	    sf_fileclose(rect[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[i],n123,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (j=0; j < n123; j++) {
		    sft[i][j] = 0;
		}
	    }

		
	    for (j=0; j < n123; j++) {
		b = ceilf(rct[i][j])+SF_ABS(sft[i][j]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }
    sf_warning("dim=%d, n123=%d",dim,n123);
	sf_warning("nn[0]=%d, nn[1]=%d, nn[2]=%d",nn[0],nn[1],nn[2]);	
	sf_warning("box[0]=%d, box[1]=%d, box[2]=%d",box[0],box[1],box[2]);
     /*added for non-stationary regularization*/
     
	/* initialize dip estimation */
	odipn_init(interp, m, n, radius, n1, n2, box, rct, sft, liter, dip0, eps, verb);

// 	odip_init(interp, m, n, radius, n1, n2, box, liter, dip0, verb);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0], n1*n2, in);
		if(slope) oslopen(wav, dip, niter, eta);
		else odipn(wav, dip, niter, eta);
		sf_floatwrite(dip[0], n1*n2, out);
	}

// 	odipn_close();
	free(dip[0]);
	free(wav[0]);
	free(dip);
	free(wav);
	return 0;
}



