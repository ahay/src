#include <stdio.h>
#include <rsf.h>

#include "chebyshev.h"

#define N 17

int main (void) 
{
    int i;
    float x, y, x1=1.0f, x2=2.0f, d[N];

    for (i=0; i < N; i++) {
	x = cosf(i*SF_PI/(N-1));
	x = 0.5*x+1.5;

	d[i] = x*x*x-x*x+1.0;
    }

    chebyshev_init(N,d,x1,x2);

    x = 1.1;
    y = x*x*x-x*x+1.0;

    printf("Compare %g and %g\n",y,chebyshev(x));
}

    
    
