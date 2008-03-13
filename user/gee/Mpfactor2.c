/* Plane prediction filter on a helix. */
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

#include <rsf.h>

#include "wilson.h"
#include "compress.h"

int main (int argc, char* argv[])
{
    int nt, nx, i, na, niter, n[3];
    float a0, s0, p, q, f, eps;
    char* lagfile;
    sf_filter ss, aa;
    sf_file filt, lag;

    sf_init (argc,argv);
    filt = sf_output("out");

    lagfile = sf_getstring("lag");
    if (NULL == lagfile) sf_error("Need lag=");
    lag = sf_output(lagfile);

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    if (!sf_getint("nx",&nx)) sf_error("Need nx=");

    n[0] = nt;
    n[1] = nx;
    n[2] = nx;

    if (!sf_getfloat("p",&p)) sf_error("Need p=");
    if (!sf_getfloat("q",&q)) sf_error("Need q=");

    ss =  sf_allocatehelix(8);
    ss->lag[0] = 1;
    ss->lag[1] = 2;
    ss->lag[2] = nt - 1;
    ss->lag[3] = nt;
    ss->lag[4] = nt + 1;
    ss->lag[5] = nt * nx - 1;
    ss->lag[6] = nt * nx;
    ss->lag[7] = nt * nx + 1;
    f = p*p*(p*p-1.) + q*q*(q*q-1.);
    s0 = 2. + 0.75*f;
    ss->flt[0] = -f;
    ss->flt[1] = 0.25*f;
    ss->flt[2] = 0.5*p*(1.-p);
    ss->flt[3] = p*p-1.;
    ss->flt[4] = -0.5*p*(1.+p);
    ss->flt[5] = 0.5*q*(1.-q);
    ss->flt[6] = q*q-1.;
    ss->flt[7] = -0.5*q*(1.+q);

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of factorization iterations */
    if (!sf_getint("na",&na)) na=25; /* filter size */
    if (25 != na) na = nt*nx + 1;
    if (!sf_getfloat("eps",&eps)) eps=FLT_EPSILON;
    /* regularization */

    aa = sf_allocatehelix (na);

    if (25 == na) {
	aa->lag[0] = 1;
	aa->lag[1] = 2;
	aa->lag[2] = nt-1;
	aa->lag[3] = nt;
	aa->lag[4] = nt+1;
	aa->lag[5] = 2*nt-1;
	aa->lag[6] = 2*nt;
	aa->lag[7] = 2*nt+1;
	aa->lag[8]  = nt*(nx-3)-2;
	aa->lag[9] = nt*(nx-3)-1;
	aa->lag[10] = nt*(nx-3);
	aa->lag[11] = nt*(nx-3)+1;
	aa->lag[12] = nt*(nx-3)+2;
	aa->lag[13] = nt*(nx-2)-2;
	aa->lag[14] = nt*(nx-2)-1;
	aa->lag[15] = nt*(nx-2);
	aa->lag[16] = nt*(nx-2)+1;
	aa->lag[17] = nt*(nx-2)+2;
	aa->lag[18] = nt*(nx-1)-1;
	aa->lag[19] = nt*(nx-1);
	aa->lag[20] = nt*(nx-1)+1;
	aa->lag[21] = nt*(nx-1)+2;
	aa->lag[22] = nt*nx-1;
	aa->lag[23] = nt*nx;
	aa->lag[24] = nt*nx+1;
    } else {
	for (i=0; i < na; i++) {
	    aa->lag[i] = i+1;
	}
    }

    for (i=0; i < na; i++) {
	aa->flt[i] = 0.;
    }

    wilson_init (nt*nx*10);
    a0 =  wilson_factor (niter, 2.*s0, ss, aa, true, 1.e-6);
    wilson_close ();

    aa = compress(aa,eps);

    if (25 != na) {
	for (i=0; i < na; i++) {
	    aa->flt[i] = 0.;
	}

	wilson_init (nt*nx*10);
	a0 =  wilson_factor (niter, 2.*s0, ss, aa, true, 1.e-6);
	wilson_close ();
    }

    na = aa->nh;

    for (i=0; i < na; i++) {
	aa->flt[i] *= a0;
    }

    sf_setformat(filt,"native_float");
    sf_putfloat(filt,"a0",a0);
    sf_putint(filt,"n1",na);
    sf_putstring(filt,"lag",lagfile);
    sf_floatwrite(aa->flt,na,filt);

    sf_setformat(lag,"native_int");
    sf_putint(lag,"n1",na);
    sf_putints(lag,"n",n,3);
    sf_intwrite(aa->lag,na,lag);

    exit(0);
}




