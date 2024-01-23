/* Shift estimation by amplitude-adjusted plane-wave destruction */
/*
  Copyright (C) 2022 University of Texas at Austin
  
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
    int n, nn, niter, order, i, liter, n1, n2, rect[3]; 
    float *u, *p, eps;
    bool verb, drift;
    sf_file inp, oth, ishift, shift;

    sf_init(argc,argv);
    inp = sf_input("in");
    oth = sf_input("other");
    shift = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float type");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) n2=1;
    n = n1*n2;
    nn = 2*n;

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* shift smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* shift smoothness on 2nd axis */
    rect[2]=1;

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */
    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    /* initialize dip estimation */
    dip3_init(n1, n2, 2, rect, liter, eps, verb);
    
    u = sf_floatalloc(nn);
    p = sf_floatalloc(nn);

    if (NULL != sf_getstring("ishift")) {
	/* initial shift */
	ishift = sf_input("ishift");
    } else {
	ishift = NULL;
    }
    
    /* read data - first monitor, then base */
    sf_floatread(u,n,oth);
    sf_floatread(u+n,n,inp);
    
    /* initialize dip */
    if (NULL != ishift) {
	sf_floatread(p,n,ishift);
	for (i=n; i < nn; i++) {
	    p[i] = 0.0f;
	}
    } else {
	for(i=0; i < nn; i++) {
	    p[i] = 0.0f;
	}	
    }
    
    /* estimate dip */
    dip3(false, 2, niter, order, 1, drift, u, p, NULL, -FLT_MAX, FLT_MAX);
	    
    /* write shift */
    sf_floatwrite(p,n,shift);

    exit(0);
}

    
