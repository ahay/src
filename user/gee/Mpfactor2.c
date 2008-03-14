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

#include "pfactor.h"

int main (int argc, char* argv[])
{
    bool fixed;
    int nt, nx, i, na, niter, n[3];
    float a0, p, q, eps;
    char* lagfile;
    sf_filter aa;
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

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of factorization iterations */
    if (!sf_getint("na",&na)) na=25; /* filter size */
    if (!sf_getfloat("eps",&eps)) eps=FLT_EPSILON;
    /* compression tolerance */
    if (!sf_getbool("fixed",&fixed)) fixed=true;
    /* if fixed size */

    if (!fixed && 25 == na) na = nt*nx + 1;

    pfactor_init(nt,nx);
    aa = sf_allocatehelix(na);

    aa = pfactor(na,p,q,niter,eps,fixed,true,&a0);
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




