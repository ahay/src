/* Stack a dataset over one of the dimensions.

Takes < gather.rsf > stack.rsf

This operation is adjoint to sfspray.
*/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int j, n1, n2, n3, i2, i3, ni, *fold = NULL, axis;
    size_t i, n;
    sf_file in, out;
    char key1[7], key2[7], *val;
    bool norm, rms;
    float *trace, *sum, f;
    sf_datatype type;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

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
    /* If y, compute the root-mean-square instead of stack */
    if (rms || !sf_getbool("norm",&norm)) norm = true;
    /* If y, normalize by fold */

    if (norm) fold = sf_intalloc (n);
    trace = sf_floatalloc (n);
    sum   = sf_floatalloc (n);
    
    for (i3=0; i3 < n3; i3++) {
	memset (sum,0,n*sizeof(float));
	if (norm) memset (fold,0,n*sizeof(int));
	
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread (trace, n, in);
	    for (i=0; i < n; i++) {
	      sum[i] += rms? trace[i]*trace[i]: trace[i];
		if (norm && (0.0 != trace[i])) fold[i]++; 
	    }
	}
	if (norm) {
	    for (i=0; i < n; i++) {
		if (fold[i] > 0) {
		  sum[i] /= fold[i];
		  if (rms) sum[i] = sqrtf(sum[i]);
		}
	    }
	}
	sf_floatwrite(sum, n, out); 
    }

    sf_close();
    exit (0);
}

/* 	$Id: stack.c,v 1.8 2004/05/06 04:43:06 fomels Exp $	 */
