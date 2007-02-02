/* Helix filters for spline in tension */
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
    const int n1=100, n[2] = {100,100};
    int i, j, na, niter;
    float s0, a0, t, eps;
    char *file;
    sf_filter ss, aa;
    sf_file flt, lag;

    sf_init(argc,argv);

    flt = sf_output("out");
    sf_setformat(flt,"native_float");

    file = sf_getstring("lag");
    if (NULL == file) sf_error("Need lag=");
    sf_putstring(flt,"lag",file);
    
    lag = sf_output(file);
    sf_setformat(lag,"native_int");    

    ss = sf_allocatehelix (12);
    ss->lag[0] = 1;
    ss->lag[1] = 2;
    i=2;
    for (j=-2; j <= 2; j++, i++) {
	ss->lag[i] = n1+j;
    }
    for (j=-2; j <= 2; j++, i++) {
	ss->lag[i] = 2*n1+j;
    }
  
    if (!sf_getfloat("tension",&t)) t=0.;
    /* spline tension */

    s0 = 36.*(114.-73.*t);
    ss->flt[0] = 96.*(8.*t-11.);
    ss->flt[1] = 84.*(1.-t);
    ss->flt[2] = 16.*(9.-8.*t);
    ss->flt[3] = 112.*(2.*t-3.);
    ss->flt[4] = ss->flt[0];
    ss->flt[5] = ss->flt[3];
    ss->flt[6] = ss->flt[2];
    ss->flt[7] = 5.*t-6.;
    ss->flt[8] = ss->flt[2];
    ss->flt[9] = ss->flt[1];
    ss->flt[10] = ss->flt[2];
    ss->flt[11] = ss->flt[7];

    aa = sf_allocatehelix (2*n1+2);
    for (i=0; i < 2*n1+2; i++) {
	aa->lag[i] = i+1;
	aa->flt[i] = 0.;
    }

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=FLT_EPSILON; 
    /* tolerance for filter compressing */

    wilson_init(20*n1);
 
    a0 = wilson_factor (niter, s0, ss, aa, true, 1.e-6);
    aa = compress (aa, eps);

    for (i=0; i < aa->nh; i++) {
	aa->flt[i] = 0.;
    }
    

    a0 = wilson_factor (niter, s0, ss, aa, true, 1.e-6);
    aa = compress (aa, eps);

    na = aa->nh;

    for (i=0; i < na; i++) {
	aa->flt[i] *= a0;
    }

    sf_putint(flt,"n1",na);
    sf_putfloat(flt,"a0",a0);

    sf_floatwrite(aa->flt,na,flt);

    sf_putint(lag,"n1",na);
    sf_putints(lag,"n",n,2);

    sf_intwrite(aa->lag,na,lag);

    exit(0);
}

