#include <rsf.h>

#include "lace.h"

#include "createhelix.h" 
#include "bound.h"
#include "pef.h"

filter lace_pef(int dim, float *dd, int jump, 
		int n, int *nd, int *center, int *gap, int *na)  
{
    int *savelags, ii[SF_MAX_DIM]; /* holding place */
    int ih, nh, lag0, j;
    filter aa;

    aa = createhelix(dim, nd, center, gap, na);  
    savelags = aa->lag;
    nh = aa->nh;

    aa->lag = sf_intalloc(nh); /* prepare interlaced helix */
    lag0 = sf_cart2line(dim, na, center);

    for (ih=0; ih < nh; ih++) {	/* sweep through the filter */
	sf_line2cart(dim, na, ih+lag0+1, ii);
	for (j=0; j < dim; j++) {
	    ii[j] -= center[j];
	}
	ii[0] *= jump;  /* interlace on 1-axis */
	aa->lag[ih] = sf_cart2line(dim, nd, ii);
    }
    na[0] *= jump;
    bound(dim, nd, nd, na, aa);  /* define  aa->mis */
    na[0] /= jump;
    
    find_pef(n, dd, aa, nh*2);    /* estimate aa coefficients */
    free(aa->lag);   
    aa->lag = savelags;		  /* restore filter lags */

    return aa;
}

