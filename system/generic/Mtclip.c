/* Clip to threshold. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2;
    float *trace, uppercut, lowercut;
    
    sf_file inp, out;

    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in inp");
    n2 = sf_leftsize(inp,1);

    /* get parameters from command line */
    if (!sf_getfloat("lowercut",&lowercut)) lowercut=0.2;
    if (!sf_getfloat("uppercut",&uppercut)) uppercut=0.8;

    trace = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	/* read data */
	sf_floatread(trace,n1,inp);

	for (i1=0; i1 < n1; i1++) {
	    if (trace[i1] < lowercut) {
		trace[i1] = 0.0;
	    } else if (trace[i1] > uppercut) {   
		trace[i1] = 1.0;
	    }
	}

	/* write result */
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}

