#include <stdio.h>

#include "kiss_fft.h"

int main(void)
{
    int n, m, m2;

    m2 = 0;
    for (n=0; n < 1001; n++) {
	m = kiss_fft_next_fast_size(n);
	if (m != m2) {
	    printf("%d ",m);
	    m2 = m;
	}
    }
    printf("\n");
    return 0;
}
