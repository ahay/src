#include <rsf.h>

#include "smoothder.h"

#include "repeat.h"
#include "causint.h"
#include "trianglen.h"
#include "weight.h"
#include "impl2.h"

static int n;
static float **tmp, *tmp1, *d;
static bool diff;

int smoothder_init(int ndim, int *rect, int *ndat, bool diff1) 
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

    tmp = sf_floatalloc2(n1,n2);
    tmp1 = tmp[0];

    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, false);    

    diff = diff1;

    if (diff) {
	impl2_init ((float) rect[0], (float) rect[1], 
		    n1, n2, 1., 50., false);
	d = sf_floatalloc(n);
    }

    return n;
}

void smoothder_close(void)
{
    free(tmp1);
    free(tmp);
    sf_conjgrad_close();
    trianglen_close();
    if (diff) {
	impl2_close();
	free(d);
    }
}

void smoothder(int niter, float* weight, float* data, float* der) 
{
    int i;

    if (NULL != weight) {
	if (diff) {
	    for (i=0; i < n; i++) {
		d[i] = data[i];
	    }
	}

	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trianglen_lop,tmp1,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trianglen_lop,tmp1,der,data,niter);
    }

    if (diff) {
	for (i=0; i < n; i++) {
	    tmp1[i] = der[i];
	}
	impl2_set(tmp);

	sf_conjgrad((NULL!=weight)? weight_lop: NULL,
		    repeat_lop,impl2_lop,tmp1,der,d,niter);
    }

    repeat_lop(false,false,n,n,der,data);
}

