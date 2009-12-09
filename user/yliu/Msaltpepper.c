/* Add salt and pepper noise to the data.*/
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

#include <math.h>
#include <time.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    float den, inten, max, *dat;
    int n1, n2, point, i, j, seed;
    bool rep, noise, allpos;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");


    if (!sf_getfloat("den",&den)) den=10.;
    /* noise density (percent, default=10, Min=0, Max=100) */
    if ( den < 0. ){
       den =0.;
    } else {
	if ( den > 100.) {
	    den = 100.;
	}
    }
    den *= 0.01;

    if (!sf_getfloat("inten",&inten)) inten=0.1;
    /* noise intensity (multiple peak value of data, default=0.1) */

    if (!sf_getbool("rep",&rep)) rep=false;
    /* if y, replace data with noise */

    if (!sf_getbool("allpos",&allpos)) allpos=false;
    /* if y, assume positive noise */

    if (!sf_getbool("noise",&noise)) noise=false;
    /* if y, output noise only */

    dat = sf_floatalloc (n1*n2);

    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */

    srand((unsigned long) seed);

    if (noise) {
        max=0.;
	sf_floatread(dat,n1*n2,in);
        for (i=0; i<n1*n2; i++) {
            if (max < fabs(dat[i])) {
               max = fabs(dat[i]);
            }
        }
        for (i=0; i<n1*n2; i++) {
             dat[i] = 0.;
        }
        for (i=0; i< n2; i++) {
            for (j=0; j<(int)(den*n1); j++) {
                point = rand()%n1;
		if (!allpos) {
		    dat[i*n1+point] = inten*max*(0.01*(rand()%100)*2.-1.);
		} else {
		    dat[i*n1+point] = inten*max*(0.01*(rand()%100));
		}
            }
        }
    } else {
        max=0.;
	sf_floatread(dat,n1*n2,in);
        for (i=0; i<n1*n2; i++) {
            if (max < fabs(dat[i])) {
               max = fabs(dat[i]);
            }
        }

        if (rep) {
            for (i=0; i< n2; i++) {
                for (j=0; j<(int)(den*n1); j++) {
                    point = rand()%n1;
		    if (!allpos) {
			dat[i*n1+point] = inten*max*(0.01*(rand()%100)*2.-1.);
		    } else {
			dat[i*n1+point] = inten*max*(0.01*(rand()%100));
		    }
                }
            }
        } else {
            for (i=0; i< n2; i++) {
                for (j=0; j<(int)(den*n1); j++) {
                    point = rand()%n1;
		    if (!allpos) {
			dat[i*n1+point] += inten*max*(0.01*(rand()%100)*2.-1.);
		    } else {
			dat[i*n1+point] += inten*max*(0.01*(rand()%100));
		    }
                }
            }
        }            
    }

    sf_floatwrite(dat,n1*n2,out);  

    exit (0);
}

/* 	$Id$	 */
