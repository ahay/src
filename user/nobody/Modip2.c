/* 2-D dip estimation by omnidirectional plane-wave destruction */
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

#include "opwd.h"

int main(int argc, char *argv[])
{
    int i, k, n1, n2, n12, n[2], niter, iter, liter, rect[2];
    float mean, a0, ai, norm, norm2, lam;
    float **dat, **ang, **p1, **p2, ***den, *dena, *rat, **out;
    sf_file inp, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n[0] = n1;
    n[1] = n2;

    if (!sf_getint("liter",&liter)) liter=100;
    /* number of linear iterations */
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* vertical smoothing */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* horizontal smoothing */

    if (!sf_getfloat("a0",&a0)) a0=0;
    /* initial angle */

    dat = sf_floatalloc2(n1,n2);
    ang = sf_floatalloc2(n1,n2);
    out = sf_floatalloc2(n1,n2);
    
    p1 = sf_floatalloc2(n1,n2);
    p2 = sf_floatalloc2(n1,n2);

    den  = sf_floatalloc3(n1,n2,2);
    dena = sf_floatalloc(n12);
    rat  = sf_floatalloc(n12);

    sf_floatread(dat[0],n12,inp);

    for (i=0; i < n12; i++) {
	ang[0][i] = a0;
	p1[0][i] = sinf(a0);
	p2[0][i] = cosf(a0);
    }

    opwd_init(n1,n2);
    sf_divn_init(2, n12, n, rect, liter, true); 

    opwd_filter(lagrange,lagrange,NULL,NULL,p1,p2,dat,out);
    
    norm = 0.;
    for (i=0; i < n12; i++) {
	out[0][i] = dat[0][i] - out[0][i];
	norm += out[0][i]*out[0][i];
    }

    for (iter=0; iter < niter; iter++) {
	sf_warning("iter=%d of %d",iter+1,niter);

	opwd_filter(lagrange_der,lagrange,NULL,NULL,p1,p2,dat,den[0]);
	opwd_filter(lagrange,lagrange_der,NULL,NULL,p1,p2,dat,den[1]);

	for(i=0; i < n12; i++) {
	    dena[i] = den[0][0][i]*p2[0][i]-den[1][0][i]*p1[0][i];
	}

        mean = 0.;
	for(i=0; i < n12; i++) {
	    mean += dena[i]*dena[i];
	}
	mean = sqrtf (n12/mean);
    
	for(i=0; i < n12; i++) {
	    out[0][i] *= mean;
	    dena[i] *= mean;
	}
    
	sf_divn (out[0],dena,rat);

	/* Choose step size */
	lam = 1.;
	for (k=0; k < 8; k++) {	    
	    for(i=0; i < n12; i++) {
		ai = ang[0][i] + lam*rat[i];
		if      (ai < -0.5*SF_PI) ai=-0.5*SF_PI;
		else if (ai >  0.5*SF_PI) ai= 0.5*SF_PI;

		p1[0][i] = sinf(ai);
		p2[0][i] = cosf(ai);
	    }
	    opwd_filter(lagrange,lagrange,NULL,NULL,p1,p2,dat,out);

	    norm2 = 0.;
	    for (i=0; i < n12; i++) {
		out[0][i] = dat[0][i] - out[0][i];
		norm2 += out[0][i]*out[0][i];
	    }
	    if (norm2 < norm) break;
	    lam *= 0.5;
	}

	for(i=0; i < n12; i++) {
	    ang[0][i] += lam*rat[i];
	    norm = norm2;
	}	
    }
    
    sf_floatwrite(ang[0],n12,dip);

    exit(0);
}
