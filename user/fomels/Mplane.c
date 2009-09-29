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
#include "contran.h"

int main(int argc, char* argv[])
{
    int i, i1, i2, ia, n1,n2, na, a1, b1, p0, nx, *lag;
    bool hyp;
    float p, s, s0, *bb, *dd, *a, *b, *pp, *qq, **aa;
    char title[10], *lagfile;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

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

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    nx = n1*n2;

    if (!sf_getbool("hyp",&hyp)) hyp=false;
    /* generate hyperbolas */

    lag = sf_intalloc(na);
    aa = sf_floatalloc2(na,nx);

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

    s = p*p; 
    s0 = p;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    bb = aa[i2*n1+i1];

	    if (hyp) {
		p = s*i2/(i1+1);
		if (p > s0) p = s0;
	    }

	    if (b1 <= 1) {
		sf_lg_int(p,na,bb);
	    } else {
		sf_taylor (p,na,dd);
		pade_apply (na,dd,a1,a,b+1);
		b[0] = 1.;
		pade_unzip(b1,b);
		pade_unzip(a1,a);
		for (i=0; i < b1-1; i++) {
		    bb[i] = b[i+1]/b[0];
		}
		for (i=b1-1, ia=0; i < na; i++, ia++) {
		    bb[i] = a[ia]/b[0];
		}
	    }
	    for (i=b1-1; i < na; i++) {
		bb[i] = - bb[i];
	    }
	}
    }

    pp = sf_floatalloc(nx);
    qq = sf_floatalloc(nx);

    contran_init (na, aa, nx, lag, true);

    sf_floatread(pp,nx,inp);
    contran_lop (false,false,nx,nx,pp,qq);
    sf_floatwrite(qq,nx,out);


    exit(0);
}



