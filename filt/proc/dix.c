#include <math.h>

#include <rsf.h>

#include "dix.h"

#include "repeat.h"
#include "causint.h"
#include "trianglen.h"
#include "weight.h"

static int n, n1, n2;
static float *tmp;

int dix_init(int ndim, int *rect, int *ndat) 
{
    int i;
    
    n=1;
    for (i=0; i <ndim; i++) {
	n *= ndat[i];
    }

    n1 = ndat[0];
    n2 = n/n1;
    
    repeat_init(n1,n2,causint_lop);
    trianglen_init(ndim,rect,ndat);

    tmp = sf_floatalloc(n);

    sf_conjgrad_init(n, n, n, n, 1., 1.e-6, true, false);    

    return n;
}

void dix_close(void)
{
    free(tmp);
    sf_conjgrad_close();
    trianglen_close();
}

void dix(int niter, float* weight, float* vrms, float* vint) 
{
    int i1, i2, i;
     
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1+i1;
	    vrms[i] *= vrms[i]*(i1+1.);
	    weight[i] /= (i1+1.); /* decrease weight with time */	    
	}
    }

    weight_init(weight);
    sf_conjgrad(weight_lop,repeat_lop,trianglen_lop,tmp,vint,vrms,niter);

    repeat_lop(false,false,n,n,vint,vrms);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1+i1;
	    vrms[i] = sqrtf(vrms[i]/(i1+1.));
	    vint[i] = sqrtf(vint[i]);
	}
    }
}
