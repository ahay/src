/* Stack a dataset over one of the dimensions.

This operation is adjoint to sfspray.
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

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int j, n1, n2, n3, i2, i3, ni, *fold = NULL, axis;
    size_t i, n;
    sf_file in, out;
    char key1[7], key2[7], *val;
    bool norm, rms, min, max;
    float *trace, *stack, f;
    double *dstack;
    sf_datatype type;

    sf_init (argc, argv);

    in  = sf_input ( "in");
    out = sf_output("out");

    if (!sf_getint("axis",&axis)) axis=2;
    /* which axis to stack */

    n1 = 1;
    for (j=0; j < axis-1; j++) {
    sprintf(key1,"n%d",j+1);
    if (!sf_histint(in,key1,&ni)) break;
    n1 *= ni;
    }

    n = (size_t) n1;

    type = sf_gettype (in);
    if (SF_FLOAT != type) {
    if (SF_COMPLEX == sf_gettype (in)) {
        n *= 2; 
	/* possibly incorrect norm for complex data */
    } else {
        sf_error("Incorrect data type in input");
    }
    }

    sprintf(key1,"n%d",axis);
    if (!sf_histint(in,key1,&n2)) sf_error("No %s= in input",key1);
    
    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
    sprintf(key1,"n%d",j+1);
    sprintf(key2,"n%d",j);
    if (!sf_histint(in,key1,&ni)) {
        sf_putint(out,key2,1);
        break;
    }
    sf_putint(out,key2,ni);
    n3 *= ni;
    
    sprintf(key1,"o%d",j+1);
    sprintf(key2,"o%d",j);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

    sprintf(key1,"d%d",j+1);
    sprintf(key2,"d%d",j);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

    sprintf(key1,"label%d",j+1);
    sprintf(key2,"label%d",j);
    if (NULL != (val = sf_histstring(in,key1))) 
        sf_putstring(out,key2,val);
    }

    if (!sf_getbool("rms",&rms)) rms = false;
    /* If y, compute the root-mean-square instead of stack. */
    if (rms || !sf_getbool("norm",&norm)) norm = true;
    /* If y, normalize by fold. */
    if (!sf_getbool("min",&min)) min = false;
    /* If y, find minimum instead of stack.  Ignores rms and norm. */
    if (!sf_getbool("max",&max)) max = false;
    /* If y, find maximum instead of stack.  Ignores rms and norm. */
    
    if (min || max) { rms = false; norm = false; }
    if (min && max) sf_error("Cannot have both min=y and max=y.");

    if (norm) fold = sf_intalloc (n);
    trace = sf_floatalloc (n);
    stack = sf_floatalloc (n);
    dstack = (min || max)? NULL: (double*) sf_alloc(n,sizeof(double));

    for (i3=0; i3 < n3; i3++) {
        if (min || max) {
            if (min) for (i=0; i<n; i++) stack[i] =  FLT_MAX;
            if (max) for (i=0; i<n; i++) stack[i] = -FLT_MAX;
        } else {
	    memset (dstack,0,n*sizeof(double));
	}
        if (norm) memset (fold,0,n*sizeof(int));
        
        for (i2=0; i2 < n2; i2++) {
            sf_floatread (trace, n, in);
            for (i=0; i < n; i++) {
                if (min || max) {
                    if (min) stack[i] = SF_MIN(stack[i],trace[i]);
                    if (max) stack[i] = SF_MAX(stack[i],trace[i]);
                } else {
		    dstack[i] += rms? (double) trace[i]*trace[i]: trace[i];
		}
                if (norm && (0.0 != trace[i])) fold[i]++; 
            }
        }
	if (!min && !max) {
	    for (i=0; i < n; i++) {
		if (norm) {
		    if (fold[i] > 0) {
			dstack[i] /= fold[i];
			if (rms) dstack[i] = sqrt(dstack[i]);
		    }
		}
		stack[i] = dstack[i];
	    }
	}
	    	
        sf_floatwrite(stack, n, out); 
    }

    exit (0);
}

