/* Inverse Q filtering by using equivalent Q value in time-frequency domain. */
/*
  Copyright (C) 2022 Jilin University
  
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
#include <math.h>

int main(int argc, char *argv[])
{
    bool  verb;
    int   n2, n1, nw, i2, i1, iw, n1w, gim;
    float d1, o1, dw, w0, w;
    float wr, det2, qq, cof;
    float *qt=NULL, *amc, *samc;
    
    kiss_fft_cpx  *phc, *opc, *cma, *outp, *asf;

    sf_file inp, out, eqt=NULL;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;
    if (!sf_getint("gim",&gim)) gim=20;    /* GIM */

    /* basic parameters */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&d1)) d1 = 1. ;
    if (!sf_histfloat(inp,"o1",&o1)) o1 = 0. ;

    n2 = sf_leftsize(inp,2);
    if(!sf_histint(inp,"n2",&nw)) sf_error("No n2=in input");
    if(!sf_histfloat(inp,"d2",&dw)) sf_error("No d2=in input");
    if(!sf_histfloat(inp,"o2",&w0)) sf_error("No o2=in input");
    sf_settype(out,SF_COMPLEX);

    if (NULL != sf_getstring("eqt")) { /* equivalent quality: eqt */
	eqt = sf_input("eqt");
	qt = sf_floatalloc(n1);
    } else {
	sf_error("Need eqt");
    }

    n1w = n1*nw;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;
    wr  = 2.*SF_PI*500;
    det2= exp(-(0.23*gim+1.63));

    asf  = (kiss_fft_cpx*)sf_complexalloc(n1w);
    cma  = (kiss_fft_cpx*)sf_complexalloc(n1w);
    outp = (kiss_fft_cpx*)sf_complexalloc(n1w);
    phc  = (kiss_fft_cpx*)sf_complexalloc(n1w);
    amc  = sf_floatalloc(n1w);
    samc = sf_floatalloc(n1w);
    opc  = (kiss_fft_cpx*)sf_complexalloc(n1w);
    
    for (i2=0; i2<n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);

	if (NULL != eqt) sf_floatread(qt,n1,eqt);

	sf_complexread((sf_complex*)asf,n1w,inp);

	for (i1=0; i1<n1; i1++) {
	    qq = qt[i1];
	    for (iw=0; iw<nw; iw++) {
		w = w0 + iw*dw;
		if (w==0) {
		    amc[i1*nw+iw] = 1.;
		    samc[i1*nw+iw] = 1.;

		    phc[i1*nw+iw].r = 1.;
		    phc[i1*nw+iw].i = 0.;
		} else {
		    cof = 1/(pow(w/wr,1/(SF_PI*qq)));
		    amc[i1*nw+iw] = exp(-1*w*i1*d1/(2*qq));
		    samc[i1*nw+iw] = (amc[i1*nw+iw]+det2)/(amc[i1*nw+iw]*amc[i1*nw+iw]+det2);

		    phc[i1*nw+iw].r = cos(cof*w*i1*d1);
		    phc[i1*nw+iw].i = sin(cof*w*i1*d1);
		}			
	    }
	}

	for (i1=0; i1<n1; i1++) {
	    for (iw=0; iw<nw; iw++) {
		opc[i1*nw+iw] = sf_cmul(asf[iw*n1+i1],phc[i1*nw+iw]);
		cma[i1*nw+iw] = sf_crmul( opc[i1*nw+iw],samc[i1*nw+iw]);
	    }
	}

	for (iw=0; iw < nw; iw++) {
	    for (i1=0; i1 < n1; i1++) {
		outp[iw*n1+i1] = cma[i1*nw+iw]; 
	    }
	}
	
	sf_complexwrite((sf_complex*)outp,n1w,out);
	
    }
}
