#include "misinput.h" 
#include "helicon.h"

void find_mask(int n, const int *known, filter aa) 
{
    int i, ih;
    float *rr, *dfre;

    rr = sf_floatalloc(n);
    dfre = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	dfre[i] = known[i]? 0.:1.;
    }
    
    helicon_init(aa);

    for (ih=0; ih < aa->nh; ih++) {
	aa->flt[ih] = 1.;
    }

    helicon_lop(false,false,n,n,dfre,rr);

    for (ih=0; ih < aa->nh; ih++) {
	aa->flt[ih] = 0.;
    }

    for (i=0; i < n; i++) {
	if ( rr[i] > 0.) aa->mis[i] = true;	
    }

    free (rr);
    free (dfre);
}
