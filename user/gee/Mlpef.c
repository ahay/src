/* Find PEF on aliased traces. */
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

#include "lace.h"
#include "printfilter.h"
#include "compress.h"

int main(int argc, char* argv[])
{
    int jump, dim, i, np;
    int n[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    float *pp;
    char varname[6], *lagfile;
    sf_filter aa;
    sf_file dat, pef, lag;

    sf_init(argc,argv);
    dat = sf_input("in");
    pef = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(pef,"lag",lagfile);

    dim = sf_filedims(dat,n);

    if (!sf_getint("jump",&jump)) jump=2;

    if (!sf_getints("a",a,dim)) sf_error("Need a=");

    if (!sf_getints("center",center,dim)) {
	for (i=0; i < dim; i++) {
	    center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
	}
    }

    if (!sf_getints("gap",gap,dim)) {
	for (i=0; i < dim; i++) {
	    gap[i] = 0;
	}
    }

    np = 1;
    for (i=0; i < dim; i++) {
	np *= n[i];
    }

    pp = sf_floatalloc(np);
    sf_floatread(pp,np,dat);

    aa = lace_pef (dim, pp, jump, np, n, center, gap, a);
    aa = compress(aa, 1.e-6);

    sf_putint(pef,"n1",aa->nh);
    sf_putint(lag,"n1",aa->nh);

    for (i=1; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	sf_putint(pef,varname,1);
	sf_putint(lag,varname,1);
	n[i] *= jump;
    }

    sf_putints (lag,"n",n,dim);
    sf_intwrite(aa->lag,aa->nh,lag);
    sf_floatwrite(aa->flt,aa->nh,pef);

    print (dim, n, center, a, aa);

    exit(0);
}

/* 	$Id$	 */
