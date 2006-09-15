#include <stdio.h>

#include "kiss_fft.h"

int main(void)
{
    int n, m, m2, m1;

    m2 = 0;
    for (n=1; n < 1001; n++) {
	m = kiss_fft_next_fast_size(n);
	if (m != m2) {
	    printf("%d ",m);
	    m2 = m;
	}
    }
    printf("\n");

    m2 = 0;
    for (n=1; n < 1001; n++) {
	m = kiss_fft_next_fast_size(n);
	if (m != m2) {
	    m1 = 2*m;
	    printf("%d ",m1);
	    m2 = m;
	}
    }
    printf("\n");

    return 0;
}
