/* Missing data interpolation in 1-D. */
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

#include "mis1.h"

int main(int argc, char* argv[])
{
    int i1, i2, n1, n2, na, niter, diter;
    float *xx, *aa;
    char *step;
    bool *known;
    sf_file in, out, filt;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

    if (NULL != sf_getstring("filtin")) {
	filt = sf_input("filtin");
	if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filtin");

	aa = sf_floatalloc(na);
	sf_floatread(aa,na,filt);
    }  else {
	filt = NULL;
	na = 2;

	aa = sf_floatalloc(2);
	aa[0] = 1.;
	aa[1] = -1.;
    }

    if (NULL == (step = sf_getstring("step"))) step = "cg";
    /* linear solver type */
    if (!sf_getint("niter",&niter)) niter=n1;
    /* number of iterations */
    if (!sf_getint("diter",&diter)) diter=niter;
    /* iteration step */

    n2 = 1+(niter-1)/diter;

    sf_putint(out,"n2",n2);
    sf_putint(out,"o2",n2==1);
    sf_putint(out,"d2",diter);

    xx = sf_floatalloc(n1);
    known = sf_boolalloc(n1);

    mis1_init(n1,na,aa);
    
    sf_floatread (xx,n1,in);
    
    for (i1=0; i1 < n1; i1++) {
	known[i1] = (bool) (xx[i1] != 0.);
    }

    for (i2=0; i2 < n2; i2++) {	
	if (1 != n2) niter = i2*diter;
        mis1 (niter, xx, known, step);
	sf_floatwrite (xx,n1,out);
    }

    exit(0);
}

/* 	$Id: Mmiss1.c 12364 2014-05-03 15:28:20Z sfomel $	 */
