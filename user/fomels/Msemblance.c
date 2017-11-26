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
    bool weighted;
    int j, n2, i2, ni, *fold = NULL, axis;
    float o2, d2;
    off_t i, n, i3, n3;
    sf_file in, out;
    char key1[7];
    float *trace;
    double t, num, den, h, sh, *dstack, *dstack2, *hstack, *hstack2;

    sf_init (argc, argv);

    in  = sf_input ( "in");
    out = sf_output("out");

    if (!sf_getbool("weighted",&weighted)) weighted=false;
    /* if use weighted semblance */

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
    sprintf(key1,"o%d",axis);
    if (!sf_histfloat(in,key1,&o2)) o2=0.0f;
    sprintf(key1,"d%d",axis);
    if (!sf_histfloat(in,key1,&d2)) d2=1.0f;
    
    n3 = sf_unshiftdim(in,out,axis);
 
    fold = sf_intalloc (n);
    trace = sf_floatalloc (n);

    dstack  = (double*) sf_alloc(n,sizeof(double));
    dstack2 = (double*) sf_alloc(n,sizeof(double));

    if (weighted) {
	hstack = (double*) sf_alloc(n,sizeof(double));
	hstack2 = (double*) sf_alloc(n,sizeof(double));
    } else {
	hstack = NULL;
	hstack2 = NULL;
    }

    for (i3=0; i3 < n3; i3++) {
	for (i=0; i < n; i++) {
	    dstack[i]  = 0.0;
	    dstack2[i] = 0.0;
	    if (weighted) {
		hstack[i]  = 0.0;
		hstack2[i] = 0.0;
	    }
	    fold[i] = 0;
	}

	sh = 0.0;
        
        for (i2=0; i2 < n2; i2++) {
	    h = o2+i2*d2;
	    sh += h;

            sf_floatread (trace, n, in);
            for (i=0; i < n; i++) {
		t = (double) trace[i];

                if (0.0 != t) {
		    dstack[i] += t;
		    dstack2[i] += t*t;

		    if (weighted) {
			hstack[i] += t*h;
			hstack2[i] += t*t*h;
		    }

		    fold[i]++; 
		}
            }
        }


	for (i=0; i < n; i++) {
	    if (fold[i] > 0) {
		if (weighted) {
		    num = (dstack[i]*hstack2[i]- dstack2[i]*hstack[i])*
			(n2*hstack[i]-dstack[i]*sh);
		    den = (n2*hstack2[i]-dstack2[i]*sh);
		    den *= den;
		} else {
		    num = dstack[i]/fold[i];
		    num *= num;
		    den = dstack2[i]/fold[i];
		}

		trace[i] = (den > 0.0f)? num/den: 0.0f;
	    } else {
		trace[i] = 0.0f;
	    }
	}
	    	
        sf_floatwrite(trace, n, out); 
    }

    exit (0);
}

