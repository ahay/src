#include <rsf.h>
#include "integrate.h"

void trapezoid(float *y, /* [n] */
	       int n,
	       float dx)
/*< 1D trapezoid >*/
/* integrate y(x) from x[0] to x[n-1] */
{
    int iy;
    float sum;

    sum=0.5f*(y[0]+y[n-1]);
    for (iy=1; iy<n-1; y++)
	sum += y[iy];
    sum *= dx;
}
