/* Display dataset signal-noise-ratio.*/
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
	int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
	int i,j,k;
        char *type;
        int ntw1, ntw2, nsw1, nsw2; /*ntw is trace-window position, nsw is sample-window position */
	float temp;  /*temporary variable*/

	float *trace;
	float *tempt; /*temporary array*/
        float en, es, snr; /*en is noise energy, es is signal energy, snr is signal-noise-ratio*/
	sf_file in;

	sf_init (argc, argv); 
	in = sf_input("in");
    
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in,2);
	/* get the trace length (n1) and the number of traces (n2) and n3*/

	if (!sf_getint("ntw1",&ntw1)) ntw1=1;
	/* trace-window beginning position (default=1)*/
	if (ntw1 < 1 || (ntw1%1)!=0.0 )  sf_error("Need positive integer input"); 
	if (!sf_getint("ntw2",&ntw2)) ntw2=n2;
	/* trace-window end position (default=n2)*/
	if (ntw2 > n2 || (ntw2%1)!=0.0 )  sf_error("Need <= n2 integer input"); 
	if (ntw1>ntw2)  sf_error("Need ntw1 <= ntw2 integer input"); 

	if (!sf_getint("nsw1",&nsw1)) nsw1=1;
	/* sample-window beginning position (default=1)*/
	if (nsw1 < 1 || (nsw1%1)!=0.0 )  sf_error("Need positive integer input"); 
	if (!sf_getint("nsw2",&nsw2)) nsw2=n1;
	/* sample-window end position (default=n1)*/
	if (nsw2 > n1 || (nsw2%1)!=0.0 )  sf_error("Need <= n1 integer input"); 
	if (nsw1>nsw2)  sf_error("Need nsw1 <= nsw2 integer input"); 

        if (NULL == (type=sf_getstring("type"))) type="stack";
        /* [stack] method type, the default is stack */

        ntw1=ntw1-1;
        ntw2=ntw2-1;
        nsw1=nsw1-1;
        nsw2=nsw2-1;

	/*set the data space*/
	trace = sf_floatalloc(n1*n2);
	tempt = sf_floatalloc(nsw2-nsw1+1);

        switch(type[0]) {
	    case 's':
                for(k=0;k<n3;k++) {
                    sf_floatread(trace,n1*n2,in);

                    en=0.0;
                    es=0.0;
                    snr=0.0;

                    for(i=0;i<nsw2-nsw1+1;i++) {
                        tempt[i]=0.0;
                    }             
                    for(i=nsw1;i<=nsw2;i++) {
                        for(j=ntw1;j<=ntw2;j++) {
                            tempt[i-nsw1]+=trace[n1*j+i];
                        }
                        es+=tempt[i-nsw1]*tempt[i-nsw1];
                    }
                    es=es/(ntw2-ntw1+1);

                    temp=0.0;
                    for(i=nsw1;i<=nsw2;i++) {
                        for(j=ntw1;j<=ntw2;j++) {
                            temp+=trace[n1*j+i]*trace[n1*j+i];
                        }
                    }
                    en=temp-es;
                    snr=10*log(es/en);
                            
                    printf("***************************************\n");
                    printf("signal energy at n3=%d      = %f \n", (k+1),es);
                    printf("noise  energy at n3=%d      = %f \n", (k+1),en);
                    printf("the SNR       at n3=%d      = %f \n", (k+1),snr);
                    printf("***************************************\n");
                }
	        break;
	    default:
	        sf_error("Unknown method type=%c",type[0]);
	        break;
        }

    exit (0);
}

/* 	$Id$	 */
