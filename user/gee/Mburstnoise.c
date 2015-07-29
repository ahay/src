/* Synthetics with bursts of noise. */
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

#include "random.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2, i, j;
    float *data, sigma, rand, thresh, thresh2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getfloat("sigma",&sigma)) sigma=1.;
    /* noise magnitude */

    if (!sf_getfloat("thresh",&thresh)) thresh=0.93;
    /* noise threshold */

    if (!sf_getfloat("thresh2",&thresh2)) thresh2=0.4;
    /* noise threshold */

    data = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	random_init (121794);
	
	for (i=0; i < n1; i++) {
	    data[i] += (random0()-0.5)*0.001;
	}
	
	random_init (121094);

	i=0;
	while (i < n1) {
	    rand = random0();
	    if (rand < thresh) { /* no noise */
		i++;
	    } else { 
		for (j=i; j < n1; j++) {
		    rand = random0();
		    data[j] += sigma*rand;
	
		    rand = random0();
		    if (rand < thresh2) break;
		}
		i = j + 1;
	    }
	}
	
	sf_floatwrite(data,n1,out);
    }


    exit(0);
}

/* 	$Id: Mburstnoise.c 4796 2009-09-29 19:39:07Z ivlad $	 */
