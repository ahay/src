#include "quantile.h"

/* changes a */
float quantile(int k, int n, float* a) {
    int i, j;
    float ak;

    ak = a[k];
    for (j=0, i=0; i < n; i++) {
	if (a[i] < ak) j++;
    }
    if (j > k) {
	for (j=0, i=0; i < n; i++) {
	    if (a[i] < ak) a[j++] = a[i];
	}
	return quantile(k, j, a);
    } else {
	for (j=0, i=0; i < n; i++) {
	    if (a[i] > ak) j++;
	}
	if (j > n - 1 - k) {
	    for (j=0, i=0; i < n; i++) {
		if (a[i] > ak) a[j++] = a[i];
	    }
	    return quantile(j + k - n, j, a);
	}
	return ak;
    }
}



