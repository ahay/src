/* 2-D multiscale dip estimation by plane wave destruction. */
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
    int n[SF_MAX_DIM], rect[3], rect2;
    int n3, i3, i1, n1, i2, n2, jt, j2, it, nt, is; 
    float p0, *u, *p, *u2, *p2, pmin, pmax, eps;
    bool verb, **mm, *m2;
    sf_file in, out, mask, idip0;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    /* get dimensions */
    dim = sf_filedims(in,n);

    n1 = n[0]; /* number of time samples */
    n2 = n[1]; /* number of traces */
    
    if (dim < 2) n2=1;

    /* make sure that n2 is a power of 2 */
    for (j=1; j < n2; j *= 2) ;
    if (j!= n2) sf_error("n2=%d, need pow(2,n)",n2);

    /* the output has twice more traces */ 
    sf_putint(out,"n2",2*n2-2);

    n123 = n1*n2;

    /* get number of 2-D frames */
    n3 = 1;
    for (j=2; j < dim; j++) {
	n3 *= n[j];
    }
    
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
    rect2 = rect[1];

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

    u  = sf_floatalloc(n123);
    p  = sf_floatalloc(n123);
    u2 = sf_floatalloc(n123);
    p2 = sf_floatalloc(n123);

    if (NULL != sf_getstring("mask")) {
	mm = sf_boolalloc2(n123,2);
	m2 = sf_boolalloc(n123);
	mask = sf_input("mask");
    } else {
	mm = (bool**) sf_alloc(2,sizeof(bool*));
	m2 = mm[0] = mm[1] = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("idip")) {
	/* initial in-line dip */
	idip0 = sf_input("idip");
    } else {
	idip0 = NULL;
    }

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("slice %d of %d;", i3+1, n3);

    	if (NULL != mask) {
	    sf_floatread(u,n123,mask);
	    mask32 (false, order, nj1, nj2, n1, n2, 1, u, mm);
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

	rect[1]=rect2;
	
	/* loop over scales */
	for (nt=n2, is=1; nt >= 2; nt /= 2, is *= 2) {

	    /* copy data */
	    for (it=i2=0; it < nt; it++, i2 += is) {
		for (i1=0; i1 < n1; i1++) {
		    jt = it*n1+i1;
		    j2 = i2*n1+i1;
		    u2[jt] = u[j2];
		    p2[jt] = is*p[j2];
		    if (NULL != m2) m2[jt] = mm[0][j2];
		}
	    }

	    /* estimate dip */
	    dip3_init(n1, nt, 1, rect, liter, eps, false);

	    /* Either scale the order or scale nj1 */

	    dip3(false, 1, niter, is*order, nj1, u2, p2, m2, is*pmin, is*pmax);
	    dip3_close();

	    /* scale dip */
	    for (i1=0; i1 < n123; i1++) p2[i1] /= is;
	    
	    /* write dip */
	    sf_floatwrite(p2,n1*nt,out);	
	    
	    /* decrease smoothness for the next scale */
	    rect[1] = SF_MAX(rect[1]/2,1);

	    /* save dip from first scale */
	    if (1==is) {
		for (i1=0; i1 < n123; i1++) p[i1] = p2[i1];
	    }
	}
    }
    sf_warning(".");
    
    exit (0);
}

/* 	$Id: Mdip2.c 7107 2011-04-10 02:04:14Z ivlad $	 */
