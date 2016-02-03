/* 2-D dip estimation by plane wave destruction. */
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

#include "dip3.h"
#include "mask6.h"

int main (int argc, char *argv[])
{
    int n123, niter, order, nj1,nj2, i,j, liter, dim;
    int n[SF_MAX_DIM], rect[3], nr, ir; 
    float p0, *u, *p, pmin, pmax, eps;
    bool verb, **mm;
    sf_file in, out, mask, idip0;

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

    n[2]= 1;
    rect[2]=1;
    nj2=1;
    
    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* dip smoothness on 2nd axis */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial dip */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("pmin",&pmin)) pmin = -FLT_MAX;
    /* minimum dip */
    if (!sf_getfloat("pmax",&pmax)) pmax = +FLT_MAX;
    /* maximum dip */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    /* initialize dip estimation */
    dip3_init(n[0], n[1], n[2], rect, liter, eps, false);

    u = sf_floatalloc(n123);
    p = sf_floatalloc(n123);

    if (NULL != sf_getstring("mask")) {
	mm = sf_boolalloc2(n123,2);
	mask = sf_input("mask");
    } else {
	mm = (bool**) sf_alloc(2,sizeof(bool*));
	mm[0] = mm[1] = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("idip")) {
	/* initial in-line dip */
	idip0 = sf_input("idip");
    } else {
	idip0 = NULL;
    }

    for (ir=0; ir < nr; ir++) {
	if (verb) sf_warning("slice %d of %d;", ir+1, nr);
    	if (NULL != mask) {
	    sf_floatread(u,n123,mask);
	    mask32 (false, order, nj1, nj2, n[0], n[1], n[2], u, mm);
	}

	/* read data */
	sf_floatread(u,n123,in);
	

	/* initialize t-x dip */
	if (NULL != idip0) {
	    sf_floatread(p,n123,idip0);
	} else {
	    for(i=0; i < n123; i++) {
		p[i] = p0;
	    }
	}
	
	/* estimate t-x dip */
	dip3(false, 1, niter, order, nj1, u, p, mm[0], pmin, pmax);
	
	/* write t-x dip */
	sf_floatwrite(p,n123,out);
	
    }
    if (verb) sf_warning(".");
    
    exit (0);
}
