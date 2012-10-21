/* Nonstationary convolution with a Ricker wavelet. Phase and Frequency can be time-varying. */
/*
  Copyright (C) 2009 China University of Petroleum-Beijing 
  and University of Texas at Austin
  
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

#include "ricker.h"

int main(int argc, char* argv[])
{
    int n1, i1, n2, i2, i, hilbn, it;
    bool norm;
    int fft_size;
    float d1, o1, mt, ee, max, freq, hilbc, *trace, *freqc, *phase, *ric, *hilb, *outtrace;
  
    sf_file in, out;

    sf_init(argc,argv);
    in  = sf_input ( "in");

    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    /* Frequency and phase*/
    if (NULL == sf_getstring("tphase") || NULL == sf_getstring("tfreq")) {
        if (!sf_getfloat("frequency",&freq)) {
			/* peak frequency for Ricker wavelet (in Hz) */
			if (!sf_getfloat("freq",&freq)) freq=0.2;
				/* peak frequency for Ricker wavelet (as fraction of Nyquist) */
		} else {
			if (!sf_histfloat(in,"d1",&d1)) d1=1.;
			freq *= 2.*d1;
        }

        trace = sf_floatalloc(n1);
        fft_size = 2*kiss_fft_next_fast_size((n1+1)/2);
        ricker_init(fft_size, 0.5*freq, 0);

        for (i2=0; i2 < n2; i2++) {
			sf_floatread(trace,n1,in);
			sf_freqfilt(n1,trace);
			sf_floatwrite(trace,n1,out);
        }
    } else { /*frequency and phase are time-varying*/
        
        sf_file tpha, tfre;
        if (!sf_histfloat(in,"o1",&o1)) o1=0;
        if (!sf_histfloat(in,"d1",&d1)) d1=0.001;
        if (!sf_getfloat("esp",&ee))  ee=0.; 
        /* if norm=y, stable parameter*/
        if (!sf_getbool("norm",&norm)) norm=false;
        tfre = sf_input ("tfreq");
        tpha = sf_input ("tphase");
        if (!sf_histint(tpha,"n1",&i) || i != n1) sf_error("the phase file has not same length");
        if (!sf_histint(tfre,"n1",&i) || i != n1) sf_error("the frequency file has not same length");
        if (!sf_getint("hiborder",&hilbn)) hilbn=6;
	/* Hilbert transformer order */
        if (!sf_getfloat("hibref",&hilbc)) hilbc=1.;
        trace = sf_floatalloc(n1);
        freqc = sf_floatalloc(n1);
        phase = sf_floatalloc(n1);
        outtrace= sf_floatalloc(n1);
        ric   = sf_floatalloc(n1);
        sf_floatread(freqc,n1,tfre);
        sf_floatread(phase,n1,tpha);

         
        sf_hilbert_init(n1, hilbn, hilbc);
        hilb = sf_floatalloc(n1);
        for (i2=0;i2 < n2; i2++) {
			sf_floatread(trace,n1,in);
			for (it=0; it < n1; it++) {
				outtrace[it] = 0.;
			}
			for (i1=0; i1 < n1; i1++) {
	
				max = 0.; 
				for (it=0;it < n1; it++) {

					mt=d1*(i1-it);
					ric[it] = (1-2*(SF_PI*freqc[i1]*mt*SF_PI*freqc[i1]*mt))*exp(-SF_PI*freqc[i1]*mt*SF_PI*freqc[i1]*mt);
	
				}
				sf_hilbert(ric,hilb);
				for (it=0; it < n1; it++) {
					hilb[it] = ric[it]*cosf(phase[i1]*SF_PI/180.) + hilb[it]*sinf(phase[i1]*SF_PI/180.);
					if (hilb[it] > max) max=hilb[it];
				}
				for (it=0;it < n1; it++) {
					if (!norm) {
						max = 1;
						ee  = 0.;
					}
					hilb[it] = hilb[it]*trace[i1]/(max+ee);
					outtrace[it] += hilb[it];
				}
			}
			sf_floatwrite(outtrace,n1,out);
        }
    }
 
    exit(0);
}
