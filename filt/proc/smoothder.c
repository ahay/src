#include <rsf.h>

#include "smoothder.h"

#include "repeat.h"
#include "causint.h"
#include "trianglen.h"
#include "weight.h"

static int n;
static float *tmp;

int smoothder_init(int ndim, int *rect, int *ndat) 
{
    int i, n1, n2;
    
    n=1;
    for (i=0; i <ndim; i++) {
	n *= ndat[i];
    }

    n1 = ndat[0];
    n2 = n/n1;
    
    repeat_init(n1,n2,causint_lop);
    trianglen_init(ndim,rect,ndat);

    tmp = sf_floatalloc(n);

    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, false);    

    return n;
}

void smoothder_close(void)
{
    free(tmp);
    sf_conjgrad_close();
    trianglen_close();
}

void smoothder(int niter, float* weight, float* data, float* der) 
{

    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trianglen_lop,tmp,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trianglen_lop,tmp,der,data,niter);
    }

    repeat_lop(false,false,n,n,der,data);
}

