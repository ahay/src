/* Prepare a filter bank for B-spline plane wave filters */
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

#include "wilson.h"
#include "compress.h"

int main(int argc, char* argv[])
{
    int nt, i, j, k, niter, np, ip, *nh, n2[2];
    float s0, p, dp, p2, pmax, p0, eps, p21;
    char *lagname, *nhname;
    sf_filter ss, aa;
    sf_file out, nhh, lag;

    sf_init (argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    nhname = sf_getstring("nh");
    if (NULL == nhname) sf_error("Need nh=");
    sf_putstring(out,"nh",nhname);

    nhh = sf_output(nhname);
    sf_setformat(nhh,"native_int");

    lagname = sf_getstring("lag");
    if (NULL == lagname) sf_error("Need lag=");
    sf_putstring(out,"lag",lagname);

    lag = sf_output(lagname);
    sf_setformat(lag,"native_int");

    if (!sf_getint("nt",&nt)) nt=40;
    /* length of the fast axis */

    ss = sf_allocatehelix (24);
    ss->lag[0] = 1; ss->lag[1] = 2; ss->lag[2] = 3;
    j = 3;
    for (k=1; k <= 3; k++) {
	for (i=-3; i <= 3; i++, j++) {
	    ss->lag[j] = k*nt+i;
	}
    }

    if (!sf_getint("np",&np)) sf_error("Need np=");
    /* number of dips */
    if (!sf_getfloat("pmax",&pmax)) pmax=2.;
    /* maximum dip */
    dp = 2.*pmax/(np-1);
    p0 = -pmax;

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=FLT_EPSILON;
    /* tolerance */

    nh = sf_intalloc(np);
    sf_putint(nhh,"n1",np);
    n2[0]=nt;
    n2[1]=20;
    sf_putints(lag,"n",n2,2);

    wilson_init (20*nt);
    for (ip=0; ip < np; ip++) {
	p = p0 + ip*dp; 
	p2 = p*p; 
	p21 = p2 + 1.;

	sf_warning("got %d of %d: %g;", ip, np, p);

	s0 = 579840.*p21;
	ss->flt[0] = 720.*(397. - 151.*p2);
	ss->flt[1] = 1152.*(25. - 151.*p2);
	ss->flt[2] =    48.*(5. - 151.*p2);
	
	ss->flt[3] =   1715.*p - 3573.*p2 - 45.;
	ss->flt[4] =  8.*(12005.*p - 10719.*p2 - 675.);
	ss->flt[5] =  5.*(84035.*p - 10719.*p21);
	ss->flt[6] = 720.*(397.*p2 - 151.);
	ss->flt[7] = -5.*(84035.*p + 10719.*p21);
	ss->flt[8] = -8.*(12005.*p + 10719.*p2 + 675.);
	ss->flt[9] = -1715.*p - 3573.*p2 - 45.;
     
	ss->flt[10] =  8.*(49.*p - 45.*p2 - 9.);
	ss->flt[11] =  64.*(343.*p - 135.*p21);
	ss->flt[12] =  8.*(12005.*p - 675.*p2 - 10719.);
	ss->flt[13] = 1152.*(25.*p2 - 151.);
	ss->flt[14] = -8.*(12005.*p + 675.*p2 + 10719.);
	ss->flt[15] = -64.*(343.*p + 135.*p21);
	ss->flt[16] = -8.*(49.*p + 45.*p2 + 9.);
     
	ss->flt[17] =  7.*p - 3.*p21;
	ss->flt[18] =  8.*(49.*p - 9.*p2 - 45.);
	ss->flt[19] =  1715.*p - 45.*p2 - 3573.;
	ss->flt[20] = 48.*(5.*p2 - 151.);
	ss->flt[21] = -1715.*p - 45.*p2 - 3573.;
	ss->flt[22] = -8.*(49.*p + 9.*p2 + 45.);
	ss->flt[23] = -7.*p - 3.*p21;

	aa = sf_allocatehelix (120);
	for (i=0; i < 120; i++) {
	    aa->flt[i] = 0.;
	}

	for (i=0; i < 25; i++) {
	    aa->lag[i] = i+1;
	}
	for (j=-30; j < 20; j++, i++) {
	    aa->lag[i] = nt+j;
	}
	for (j=-24; j < 6; j++, i++) {
	    aa->lag[i] = 2*nt+j;
	}
	for (j=-11; j < 4; j++, i++) {
	    aa->lag[i] = 3*nt+j;
	}

	(void) wilson_factor (niter, s0, ss, aa, false, 1.e-6);
	aa = compress (aa, eps);

	for (i=0; i < aa->nh; i++) {
	    aa->flt[i] = 0.;
	}

	(void) wilson_factor (niter, s0, ss, aa, false, 1.e-6);

	sf_floatwrite (aa->flt,aa->nh,out);
	sf_intwrite (aa->lag,aa->nh,lag);

	nh[ip] = aa->nh;
	sf_deallocatehelix (aa);
    }
    sf_warning(".");

    sf_intwrite (nh,np,nhh);

    exit(0);
}
