/* Compute data envelope.
*/
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

#include "triangle.h"

int main (int argc, char* argv[])
{
    bool freq;
    int n1,n2,n3, i1,i2,i3, tc1, tc2;
    float *data, **bot=NULL, **top=NULL, den;
    float complex *cdat, *ctop=NULL;
    triangle tr1=NULL, tr2=NULL;
    kiss_fft_cfg forw, invs;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1; 
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("freq",&freq)) freq=false;
    /* if y, compute instantenous frequency */

    data = sf_floatalloc(n1);
    cdat = sf_complexalloc (n1);

    if (freq) {
	ctop = sf_complexalloc (n1);
	top = sf_floatalloc2(n1,n2);
	bot = sf_floatalloc2(n1,n2);

	if (!sf_getint("tc1",&tc1)) tc1=1;
	if (!sf_getint("tc2",&tc2)) tc2=1;
	/* smoothing triangle size for instanteneous frequency (if freq=y) */
    
	tr1 = triangle_init (tc1,n1);
	tr2 = triangle_init (tc2,n2);
    }

    forw = kiss_fft_alloc(n1,0,NULL,NULL);
    invs = kiss_fft_alloc(n1,1,NULL,NULL);
    if (NULL == forw || NULL == invs) sf_error("KISS FFT allocation error");

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(data,n1,in);
	    for (i1=0; i1 < n1; i1++) {
		cdat[i1]=data[i1];
	    }
	    kiss_fft(forw,(const kiss_fft_cpx *) cdat, 
		     (kiss_fft_cpx *) cdat);
	    cdat[0] *= 0.5;
	    cdat[n1/2] *= 0.5;
	    for (i1=n1/2+1; i1 < n1; i1++) {
		cdat[i1]=0.;
	    }
	    if (freq) {
		for (i1=0; i1 < n1; i1++) {
		    ctop[i1] = (2.*SF_PI*I*i1/n1) * cdat[i1];
		}
		kiss_fft(invs,(const kiss_fft_cpx *) ctop,
			 (kiss_fft_cpx *) ctop);
		kiss_fft(invs,(const kiss_fft_cpx *) cdat,
			 (kiss_fft_cpx *) cdat);
		for (i1=0; i1 < n1; i1++) {
		    bot[i2][i1] = crealf(conjf(cdat[i1])*cdat[i1]);
		    top[i2][i1] = crealf(-conjf(cdat[i1])*ctop[i1]*I);
		}
		smooth(tr1,i2*n1,1,false,top[0]);
		smooth(tr1,i2*n1,1,false,bot[0]);
	    } else {
		kiss_fft(invs,(const kiss_fft_cpx *) cdat,
			 (kiss_fft_cpx *) cdat);
		for (i1=0; i1 < n1; i1++) {
		    data[i1] = 2.*cabsf(cdat[i1])/n1;
		}
		sf_floatwrite(data,n1,out);
	    }
	}
	if (freq) {
	    smooth(tr2,0,n1,false,bot[0]);
	    smooth(tr2,0,n1,false,top[0]);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    den = bot[i2][i1];
		    if (den != 0.) top[i2][i1] /= den;
		}
	    }
	    sf_floatwrite(top[0],n1*n2,out);
	}
    }

    exit(0);
}

/* 	$Id$	 */
