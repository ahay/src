#include "quantile.h"

/* changes a */
float quantile(int q, int n, float* a) {
    float *i, *j, ak, *low, *hi, buf, *k;

    low=a;
    hi=a+n-1;
    k=a+q; 
    while (low<hi) {
	ak = *k;
	i = low; j = hi;
	do {
	    while (*i < ak) i++;     
	    while (*j > ak) j--;     
	    if (i<=j) {
		buf = *i;
		*i++ = *j;
		*j-- = buf;
	    }
	} while (i<=j);
	if (j<k) low = i; 
	if (k<i) hi = j;
    }
    return (*k);
}
