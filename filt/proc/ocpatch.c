#include <stdio.h>

#include <rsf.h>

#include "ocpatch.h"

static int n1, n2;
static long **table; 

void ocpatch_init(int dim, int nw, int np, int* npatch, int* nwall, int* nwind)
{
    int i, ip, i2, ff[SF_MAX_DIM], t2, t;
    
    n1 = nwind[0];
    n2 = nw/n1;

    table = (long**) sf_alloc(np,sizeof(long*));
    table[0] = (long*) sf_alloc(np*n2,sizeof(long));

    for (ip=0; ip < np; ip++) {
	if (ip) table[ip] = table[0]+ip*n2;

	t = ip;
	for(i = 0; i < dim; i++) {
	    t2 = t%npatch[i];
	    t /= npatch[i];

	    if(npatch[i] == 1) {
		ff[i] = 0;
	    } else if (t2 == npatch[i]-1) {
		ff[i] = nwall[i] - nwind[i];
	    } else {	    
		ff[i] = t2*(nwall[i] - nwind[i])/(npatch[i] - 1.0);
	    }
	}	
	t = sf_cart2line (dim-1, nwall+1, ff+1);
	table[ip][0] = (t*nwall[0] + ff[0])*sizeof(float);
	for (i2=1; i2 < n2; i2++) {
	    t = i2;
	    for (i = 1; i < dim-1; i++) {
		/* cartesian coordinates in window */
		ff[i] += t%nwind[i];
		t /= nwind[i];
	    }
	    t += ff[dim-1];
	    for (i = dim-2; i >= 1; i--) {
		/* line coordinates in input */
		t = t*nwall[i] + ff[i];
	    }
	    table[ip][i2] = (t*nwall[0] + ff[0])*sizeof(float); 
	}
    }
}

void ocpatch_close(void)
{
    free(table[0]);
    free(table);
}

void ocpatch_lop (int ip, bool adj, FILE *wall, float* wind)
{
    int i2;

    for (i2=0; i2 < n2; i2++, wind += n1) {
	if (0 != fseek(wall,table[ip][i2],SEEK_SET))
	    sf_error("%s: seeking error:",__FILE__);

	if (adj) {
	    if (n1 != fwrite(wind,sizeof(float),n1,wall))
		sf_error("%s: writing error:",__FILE__);
	} else {
	    if (n1 != fread(wind,sizeof(float),n1,wall))
		sf_error("%s: reading error:",__FILE__);
	}
    }
}
