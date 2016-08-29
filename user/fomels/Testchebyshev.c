#include <stdio.h>
#include <time.h>
#include <rsf.h>

#include "chebyshev.h"

#define N 9

int main (void) 
{
    int i, seed;
    float x, y, x1=1.0f, x2=2.0f, d[N];

    seed = time(NULL);
    init_genrand((unsigned long) seed);

    chebyshev_init(N,x1,x2);

    for (i=0; i < N; i++) {
	x = cosf(i*SF_PI/(N-1));
	x = 0.5*x+1.5;

	d[i] = (x*x*x-x*x+1.0);
    }

    chebyshev_set(d);

    x = 1.0+genrand_real1();
    y = (x*x*x-x*x+1.0);

    printf("Compare %g and %g\n",y,chebyshev(x));

    chebyshev_close();

    return 0;
}

    
    
