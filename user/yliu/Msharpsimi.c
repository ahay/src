/* Sharpen similarity measure between two datasets. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
    int n1, n2, n3, n12, sn1, sn2, sn3, i1, i2, nit, niter; 
    bool verb;    
    float perc, eps, *t1, *t2, *xk, * yk, *simi;
    sf_file in, out, other;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    other = sf_input("other");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histint(other,"n1",&sn1)) sf_error("No n1= in other");
    if (!sf_histint(other,"n2",&sn2)) sf_error("No n2= in other");

    n3 = sf_leftsize(in,2);
    sn3 = sf_leftsize(other,2);

    if (n1!=sn1 || n2!=sn2 || n3!=sn3) sf_error("Different size in two datasets");

    n12 = n1*n2;
    n3 = sf_leftsize(in,2);
    sn3 = sf_leftsize(other,2);
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("perc",&perc)) perc=98.;
    /* percentage for sharpen, default is 98*/

    if (perc < 0. || perc > 100.)  sf_error("Need input is in range [0.0,100.0]"); 

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    t1 = sf_floatalloc(n12);
    t2 = sf_floatalloc(n12);
    xk = sf_floatalloc(n12);
    yk = sf_floatalloc(n12);
    simi = sf_floatalloc(n12);

    sf_sharpen_init(n12,perc,0.5);

    for (i2=0; i2 < n3; i2++) {
	if (verb) sf_warning("Slice %d of %d",i2+1,n2);
	sf_floatread(t1,n12,in);
	sf_floatread(t2,n12,other);

	for(i1=0; i1 < n12; i1++) {
	    xk[i1] = 0.;
	    yk[i1] = 0.;
	}
	for(nit=0; nit < niter; nit++) {
	    if (verb) sf_warning("Iteration %d of %d",nit+1,niter);
	    for(i1=0; i1 < n12; i1++) {
		xk[i1] = ((t1[i1]-t2[i1])*t2[i1] + eps*eps*xk[i1])/(t2[i1]*t2[i1]+eps*eps);
		yk[i1] = ((t2[i1]-t1[i1])*t1[i1] + eps*eps*yk[i1])/(t1[i1]*t1[i1]+eps*eps);
/*		xk[i1] = (t1[i1]*t2[i1] + eps*eps*xk[i1])/(t2[i1]*t2[i1]+eps*eps);
		yk[i1] = (t2[i1]*t1[i1] + eps*eps*yk[i1])/(t1[i1]*t1[i1]+eps*eps);
 */	    }
	    sf_sharpen(xk);
	    sf_weight_apply(n12,xk);
	    sf_sharpen(yk);
	    sf_weight_apply(n12,yk);
	}
	for(i1=0; i1 < n12; i1++) {
	    simi[i1] = (1+xk[i1])*(1+yk[i1]);
/*	    simi[i1] = xk[i1]*yk[i1];
 */	}
	sf_floatwrite (simi,n12,out);	
    }

    exit (0);
}

/* 	$Id$	 */

