/* Generating plane waves with steering filters. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <stdio.h>

#include <rsf.h>

#include "pade.h"

int main(int argc, char* argv[])
{
    int i, i1, i2, ia, n1,n2, na, a1, b1, p0, nx, *lag, n[2];
    bool hyp;
    float p, s, s0, *bb, *dd, *a, *b;
    char title[10], *lagfile;
    sf_file in, out, lg;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    lg = sf_output("lag");

    if (!sf_getfloat("p",&p)) p=0.7;
    /* plane wave slope */

    if (!sf_getint("a1",&a1)) a1=2;
    /* filter length */

    if (!sf_getint("b1",&b1)) b1=1;
    /* denominator length */

    snprintf(title,9,"%d points",a1);
    sf_putstring(out,"title",title);

    if (NULL != (lagfile = sf_getstring("lag"))) 
	sf_putstring(out,"lag",lagfile);

    na = a1+b1-1;
    p0 = (a1-1)/2;

    if (!sf_getint("n1",&n1)) sf_error("No n1= in input");
    if (!sf_getint("n2",&n2)) sf_error("No n2= in input");

    if (!sf_getbool("hyp",&hyp)) hyp=false;
    /* generate hyperbolas */

    nx = n1*n2;
    lag = sf_intalloc(na);
    bb = sf_floatalloc(na);

    if (b1 > 1) {
	a = sf_floatalloc(a1);
	b = sf_floatalloc(b1);
	dd = sf_floatalloc(na);
	pade_init(b1-1);
    } else {
	a = NULL;
	b = NULL;
	dd = NULL;
    }

    for (ia=0; ia < b1-1; ia++) {
	lag[ia] = ia+1;
    }
    for (ia=b1-1; ia < na; ia++) {
	lag[ia] = n1+ia+1-p0-b1;
    }
    n[0] = na;
    n[1] = nx;

    sf_setformat(lg,"native_int");
    sf_putints(lg,"n",n,2);
    sf_putint(lg,"n1",na);
    sf_putint(lg,"n2",1);
    sf_intwrite(lag,na,lg);

    sf_setformat(out,"native_float");
    sf_putint(out,"n1",na);
    sf_putint(out,"n2",nx);

    s = p*p; 
    s0 = p;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (hyp) {
		p = s*i2/(i1+1);
		if (p > s0) p = s0;
	    }

	    if (b1 <= 1) {
		sf_lg_int(p,na,bb);
	    } else {
		sf_taylor (p,na,dd);
		pade_apply (na,dd,a,b+1);
		b[0] = 1.;
		pade_unzip(b1,b);
		pade_unzip(a1,a);
		for (i=0; i < b1-1; i++) {
		    bb[i] = b[i+1]/b[0];
		}
		for (i=b1-1, ia=0; i < na; i++, ia++) {
		    b[i] = a[ia]/b[0];
		}
	    }
	    for (i=b1-1; i < na; i++) {
		bb[i] = - bb[i];
	    }
	    
	    sf_floatwrite(bb,na,out);
	}
    }
    
    exit(0);
}



