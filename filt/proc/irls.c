#include <math.h>

#include <rsf.h>

#include "irls.h"

#include "quantile.h"

static float *abs1;

void irls_init(int n) {
    abs1 = sf_floatalloc(n);
}

void irls_close(void) {
    free (abs1);
}

void l1 (int n, const float *res, float *weight)  {
    float rbar;
    int i;

    for (i=0; i < n; i++) {
	abs1[i] = fabsf(res[i]);
    }

    rbar = quantile(n/2,n,abs1);

    sf_warning("in l1");

    for (i=0; i < n; i++) {
	weight[i] = 1./sqrtf(hypotf(1.,res[i]/rbar));
    }
}

void cauchy (int n, const float *res, float *weight)  {
    float rbar;
    int i;

    for (i=0; i < n; i++) {
	abs1[i] = fabsf(res[i]);
    }

    rbar = quantile(n/2,n,abs1);

    sf_warning("in cauchy");

    for (i=0; i < n; i++) {
	weight[i] = 1./hypotf(1.,res[i]/rbar);
    }
}
