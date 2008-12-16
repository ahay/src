#include <stdlib.h>
#include <stdio.h>

#include "quantile.h"

#define NX 5

int main (int argc, char* argv[])
{
    int i, j;
    float xx[] = {1.,3.,4.,5.,2.};
    float yy[NX];
    
    for (i=0; i < NX; i++) {
	for (j=0; j < NX; j++) {
	    yy[j] = xx[j];
	}
	printf("%d %f\n", i, sf_quantile (i, NX, yy));
    }
    
    exit (0);
}
