/* B-spline plane-wave filter */
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

int main(int argc, char* argv[])
{
    const int n1=100, n2[2]={100,100};
    int i, niter, nw, na;
    float s0, a0, p, p2, p21, eps;
    char *lagfile;
    sf_filter ss, aa;
    sf_file out, lag;

    sf_init(argc,argv);
    out = sf_output("out");
    lag = sf_output("lag");
    sf_setformat(out,"native_float");
    sf_setformat(lag,"native_int");

    if (NULL != (lagfile = sf_getstring("lag"))) 
	sf_putstring(out,"lag",lagfile);

    if (!sf_getint("nw",&nw)) nw=2; /* filter size */
    if (!sf_getfloat("p",&p)) p=0.; /* plane-wave slope */

    p2 = p*p; 
    p21 = p2 + 1.;

    switch(nw) {
	case 2:
	    ss = sf_allocatehelix(4);
	    ss->lag[0]=1;
	    ss->lag[1]=n1-1;
	    ss->lag[2]=n1;
	    ss->lag[3]=n1+1;

	    s0 = 8.*p21;
	    ss->flt[0] = 2.*(1.-2.*p2);
	    ss->flt[1] =  3.*p-p21;
	    ss->flt[2] = 2.*(p2-2.);
	    ss->flt[3] = -3.*p-p21;
	    break;
	case 3:
	    ss = sf_allocatehelix(12);
	    ss->lag[0]=1;
	    ss->lag[1]=2;
	    ss->lag[2]=n1-2;
	    ss->lag[3]=n1-1;
	    ss->lag[4]=n1;
	    ss->lag[5]=n1+1;
	    ss->lag[6]=n1+2;
	    ss->lag[7]=2*n1-2;
	    ss->lag[8]=2*n1-1;
	    ss->lag[9]=2*n1;
	    ss->lag[10]=2*n1+1;
	    ss->lag[11]=2*n1+2;

	    s0 = 792.*p21;
	    ss->flt[0] = 24.*(13. - 11.*p2);
	    ss->flt[1] = 12.*( 1. - 11.*p2);

	    ss->flt[2] =  2.*(25.*p - 2. - 26.*p2);
	    ss->flt[3] =  4.*(125.*p - 26.*p21);
	    ss->flt[4] = 24.*(13.*p2 - 11.);
	    ss->flt[5] = -4.*(125.*p + 26.*p21);
	    ss->flt[6] = -2.*(25.*p + 2. + 26.*p2);

	    ss->flt[7]  =  5.*p - 2.*p21;
	    ss->flt[8]  =  2.*(25.*p - 26. - 2.*p2);
	    ss->flt[9] = 12.*(p2 - 11.);
	    ss->flt[10] = -2.*(25.*p + 26. + 2.*p2);
	    ss->flt[11] = -5.*p - 2.*p21;
	    break;
	case 4:
	    ss = sf_allocatehelix(24);

	    ss->lag[0]=1;
	    ss->lag[1]=2;
	    ss->lag[2]=3;
	    ss->lag[3]=n1-3;
	    ss->lag[4]=n1-2;
	    ss->lag[5]=n1-1;
	    ss->lag[6]=n1;
	    ss->lag[7]=n1+1;
	    ss->lag[8]=n1+2;
	    ss->lag[9]=n1+3;
	    ss->lag[10]=2*n1-3;
	    ss->lag[11]=2*n1-2;
	    ss->lag[12]=2*n1-1;
	    ss->lag[13]=2*n1;
	    ss->lag[14]=2*n1+1;
	    ss->lag[15]=2*n1+2;
	    ss->lag[16]=2*n1+3;
	    ss->lag[17]=3*n1-3;
	    ss->lag[18]=3*n1-2;
	    ss->lag[19]=3*n1-1;
	    ss->lag[20]=3*n1;
	    ss->lag[21]=3*n1+1;
	    ss->lag[22]=3*n1+2;
	    ss->lag[23]=3*n1+3;

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
	    break;
	default:
	    ss = NULL;
	    s0 = 0.;
	    sf_error("No support for size nw=%d",nw);
	    break;
    }

    na = (nw-1)*(n1+1);
    aa = sf_allocatehelix(na);
    for (i=0; i < na; i++) {
	aa->lag[i] = i+1;
	aa->flt[i] = 0.0;
    }

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of spectral decomposition iterations */
    if (!sf_getfloat("eps",&eps)) eps=SF_EPS;

    wilson_init(20*n1);
    a0 = wilson_factor (niter, s0, ss, aa, true, 1.e-6);	
    aa = compress(aa, eps);
  
    for (i=0; i < aa->nh; i++) {
	aa->flt[i] = 0.;
    }

    a0 = wilson_factor (niter, s0, ss, aa, true, 1.e-6);	
    aa = compress(aa, eps);
    na = aa->nh;

    for( i=0; i < na; i++) {
	aa->flt[i] *= a0;
    }

    sf_putint(out,"n1",na);
    sf_putfloat(out,"a0",a0);
    sf_floatwrite(aa->flt,na,out);

    sf_putint(lag,"n1",na);
    sf_putints(lag,"n",n2,2);
    sf_intwrite(aa->lag,na,lag);

    exit(0);
}

