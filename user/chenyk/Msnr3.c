/* Compute signal-noise-ratio.
SNR=10 log10(sum(clean)/sum(noise))*/
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
	int n1,nn1,n2,nn2; /*n1 is trace length, n2 is the number of traces*/
	int i,j;

	float *trace1, *trace2;
        float en, es, snr; /*en is noise energy, es is signal energy, snr is signal-noise-ratio*/
	sf_file signal, noise, snrf;

	sf_init (argc, argv); 
	signal = sf_input("in");
	noise = sf_input("noise");
        snrf = sf_output("out");
    
	if (!sf_histint(signal,"n1",&n1)) sf_error("No n1= in input");
	n2 = sf_leftsize(signal,1);
	/* get the trace length (n1) and the number of traces (n2) from signal*/

	if (!sf_histint(noise,"n1",&nn1)) sf_error("No n1= in noise");
	nn2 = sf_leftsize(noise,1);
	/* get the trace length (n1) and the number of traces (n2) from noise*/	

	if(n1 !=nn1 || n2 !=nn2) sf_error("size doesn't match");

        sf_putint(snrf,"n1",1);
	sf_putfloat(snrf,"d1",1);
        sf_putint(snrf,"n2",1);
        sf_putint(snrf,"n3",1);

	/*set the data space*/
	trace1 = sf_floatalloc(n1*n2);
	trace2 = sf_floatalloc(n1*n2);

                    sf_floatread(trace1,n1*n2,signal);
		    sf_floatread(trace2,n1*n2,noise);

                    en=0.0;
                    es=0.0;
                    snr=0.0;
       
                    for(i=0;i<=n1-1;i++) {
                        for(j=0;j<=n2-1;j++) {
                            es+=trace1[n1*j+i]*trace1[n1*j+i];
                        }
                    }

                    for(i=0;i<=n1-1;i++) {
                        for(j=0;j<=n2-1;j++) {
                            en+=trace2[n1*j+i]*trace2[n1*j+i];
                        }
                    }
                    snr=10*log10(es/en);

		    sf_floatwrite(&snr,1,snrf);

    exit (0);
}
