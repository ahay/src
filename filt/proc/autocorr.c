
#include "autocorr.h"
#include "helix.h"
#include "compress.h"

filter autocorr(const filter aa, float a0, float *s0, float eps)
{
    int i, j, k, n, na;
    float f, b0;
    filter ss;

    na = aa->nh;

    ss = allocatehelix (na*(na+1)/2);
 
    b0 = a0*a0;
    for (i=0; i < na; i++) {
	f = aa->flt[i];
	b0 += f*f;
	ss->flt[i] = a0*f;
	ss->lag[i] = aa->lag[i];
    }
    *s0 = b0;
	  
    k = na-1;
    for (i=0; i < na; i++) {
	for (j = i+1; j < na; j++) {
	    k++;

	    ss->flt[k] = aa->flt[j] * aa->flt[i];
	    ss->lag[k] = aa->lag[j] - aa->lag[i];

	    for (n=0; n < k-1; n++) {
		if (ss->lag[n] == ss->lag[k] && ss->flt[n] != 0.) {
		    ss->flt[n] += ss->flt[k];
		    ss->flt[k] = 0.;
		}
	    }
	}
    }

    return compress(ss,eps);
}
