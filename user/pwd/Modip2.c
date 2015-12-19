/* 2-D dip estimation by omnidirectional plane wave destruction. */
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

#include "odip2.h"

int main (int argc, char *argv[])
{
    int n123, niter, order, i,j, liter, dim;
    int n[SF_MAX_DIM], rect[2], nr, ir; 
    float a0, *u, *a, eps;
    bool verb;
    sf_file in, out, angle;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    dim = sf_filedims(in,n);
    if (dim < 2) n[1]=1;
    n123 = n[0]*n[1];
    nr = 1;
    for (j=2; j < dim; j++) {
	nr *= n[j];
    }
    
    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* dip smoothness on 2nd axis */

    if (!sf_getfloat("a0",&a0)) a0=0.;
    /* initial dip */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization */

    /* initialize dip estimation */
    odip2_init(n[0], n[1], rect, liter, verb);

    u = sf_floatalloc(n123);
    a = sf_floatalloc(n123);

    if (NULL != sf_getstring("angle")) {
	/* initial in-line dip */
	angle = sf_input("angle");
    } else {
	angle = NULL;
    }

    for (ir=0; ir < nr; ir++) {

	/* read data */
	sf_floatread(u,n123,in);
	

	/* initialize t-x dip */
	if (NULL != angle) {
	    sf_floatread(a,n123,angle);
	} else {
	    for(i=0; i < n123; i++) {
		a[i] = a0;
	    }
	}
	
	/* estimate dip */
	odip2(niter, order, u, a, eps);
	
	/* write dip */
	sf_floatwrite(a,n123,out);
    }
    
    exit (0);
}

