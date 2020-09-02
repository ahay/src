#include <assert.h> /* for assert */
#include <float.h>  /* for DBL_EPSILON */

int main(void)
{
    int i;
    double eps, one;

    eps = 1.0;
    for (i=0; i < 100; i++) {
	eps /= 2;
	one = 1.0+eps;

	/* !!! INSERT SOMETHING HERE !!! */
    }

    assert(DBL_EPSILON==eps);
}
