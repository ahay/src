#include <stdlib.h>
#include <alloca.h>
#include <stdio.h>

#include <rsf.h>

#include "printfilter.h"
#include "boxfilter.h"
#include "helix.h"

void print (int dim, const int *nd, const int *center, const int *na, 
	    const filter aa) 
{
    float* filt; 
    int ii[SF_MAX_DIM], ii0[SF_MAX_DIM], ia, i, j, k, nf, br; 

    nf = 1;
    for (j=0; j < dim; j++) {
	nf *= na[j];
    }
    filt = sf_floatalloc(nf);

    box (dim, nd, center, na, aa, nf, filt); 
    fprintf(stderr,"-----------------------------------------\n");

    for (j=0; j < dim; j++) {
	ii0[j] = 0;
    }

    for (ia = 0; ia < nf - na[0]+1; ia += na[0]) {
	sf_line2cart(dim, na, ia, ii);
	br = 0;
	for (j=0; j < dim; j++) {
	    if (ii[j] != ii0[j]) br++;
	    ii0[j] = ii[j];
	}

        for (i=0; i < br-1; i++)
	    fprintf (stderr,"+++++++++++++++++++++++++++++++++++++++++\n");

	for (k=0; k < na[0]; k++) {
	    fprintf(stderr,"%7.3f",filt[k+ia]);
	}
	fprintf (stderr,"\n");
    }
    
    fprintf (stderr,"-----------------------------------------\n");

    free (filt);
}

