#include <rsf.h>

#include "smoothder.h"

#include "repeat.h"
#include "causint.h"
#include "trianglen.h"
#include "trisl.h"
#include "weight.h"
#include "impl2.h"

static int n, n2;
static float **tmp, *tmp1, *tmp2;
static bool diff, dip;

int smoothder_init(int ndim, int *rect, int *ndat, bool diff1, bool dip1) 
{
    int i, n1;
    
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
    tmp2 = sf_floatalloc(n);

    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, diff);    

    diff = diff1;
    dip = dip1;

    if (dip) {
	trisl_init(n1,n2,rect[0],rect[1]);
    } else if (diff) {
	impl2_init ((float) rect[0], (float) rect[1], 
		    n1, n2, 1., 50., false);
    }

    return n;
}

void smoothder_close(void)
{
    free(tmp1);
    free(tmp);
    sf_conjgrad_close();
    trianglen_close();
    if (diff) impl2_close();
    if (dip) trisl_close();
}

void smoothder(int niter, float* weight, float* data, float* der) 
{
 
    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trianglen_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trianglen_lop,tmp2,der,data,niter);
    }

    repeat_lop(false,false,n,n,der,data);
}

void smoothdip(int niter, float** dip, float* weight, float* data, float* der) 
{
    trisl_set(dip);

    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trisl_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trisl_lop,tmp2,der,data,niter);
    }
    
    repeat_lop(false,false,n,n,der,data);
}


void smoothdiff(int niter, int ncycle, float* weight, float* data, float* der) 
{
    int i, iter;

    for (i=0; i < n; i++) {
	tmp2[i] = 0.;
    }

    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,impl2_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,impl2_lop,tmp2,der,data,niter);
    }
	
    for (iter=0; iter < ncycle; iter++) {
	/*
	for (i=0; i < n; i++) {
	    tmp1[i] = der[i];
	}
	impl2_set(tmp);
	*/
	sf_conjgrad((NULL!=weight)? weight_lop: NULL,
		    repeat_lop,impl2_lop,tmp2,der,data,niter/ncycle);
    }    

    repeat_lop(false,false,n,n,der,data);
}

