#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "grid1.h"

int main(void) {
    int i;
    float x[] = {1., 5., 2, 3., 10.};
    float z[] = {0., 1., 2., 5., 4.5};
    float y[2];
    grid1 grid;

    grid = grid1_init();

    for (i=0; i < 5; i++) {
		y[0] = sinf(x[i]);
		y[1] = cosf(x[i]);
		grid1_insert(grid,x[i],2,y);
    }

    for (i=0; i < 5; i++) {
		grid1_interpolate (grid,x[i],2,y);
		printf("%g = %g\n%g = %g\n\n",y[0],sinf(x[i]),y[1],cosf(x[i]));

		grid1_interpolate (grid,z[i],2,y);
		printf("sin(%g) = %g\ncos(%g) = %g\n\n",z[i],y[0],z[i],y[1]);
    }
    
    grid1_close(grid);

    exit(0);
}
