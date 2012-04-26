/* Semblance over the specified axis. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
    int j, n2, i2, ni, *fold = NULL, axis;
    off_t i, n, i3, n3;
    sf_file in, out;
    char key1[7];
    float *trace;
    double t, *dstack, *dstack2;

    sf_init (argc, argv);

    in  = sf_input ( "in");
    out = sf_output("out");

    if (!sf_getint("axis",&axis)) axis=2;
    /* which axis to stack */

    n = 1;
    for (j=0; j < axis-1; j++) {
	sprintf(key1,"n%d",j+1);
	if (!sf_histint(in,key1,&ni)) break;
	n *= ni;
    }

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    sprintf(key1,"n%d",axis);
    if (!sf_histint(in,key1,&n2)) sf_error("No %s= in input",key1);
    n3 = sf_unshiftdim(in,out,axis);
 
    fold = sf_intalloc (n);
    trace = sf_floatalloc (n);

    dstack  = (double*) sf_alloc(n,sizeof(double));
    dstack2 = (double*) sf_alloc(n,sizeof(double));

    for (i3=0; i3 < n3; i3++) {
	for (i=0; i < n; i++) {
	    dstack[i]  = 0.0;
	    dstack2[i] = 0.0;
	    fold[i] = 0;
	}
        
        for (i2=0; i2 < n2; i2++) {
            sf_floatread (trace, n, in);
            for (i=0; i < n; i++) {
		t = (double) trace[i];

                if (0.0 != t) {
		    dstack[i] += t;
		    dstack2[i] += t*t;
		    fold[i]++; 
		}
            }
        }


	for (i=0; i < n; i++) {
	    if (fold[i] > 0) {
		dstack[i]  /= fold[i];
		dstack2[i] /= fold[i];
	    }
	    trace[i] = dstack[i]*dstack[i]/dstack2[i];
	}
	    	
        sf_floatwrite(trace, n, out); 
    }

    exit (0);
}

