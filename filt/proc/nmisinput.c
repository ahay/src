#include "nmisinput.h"
#include "nhelicon.h"


void find_mask(int nd, const int *known, nfilter aa) {
    float *rr, *dfre;
    int ip, i;

    rr = sf_floatalloc(nd);
    dfre = sf_floatalloc(nd);

    for (i=0; i < nd; i++) {
	dfre[i] = known[i]? 0.:1.;
    }

    nhelicon_init( aa);

    for (ip=0; ip < aa->np; ip++) {
	for (i=0; i < aa->hlx[ip]->nh; i++) {
	    aa->hlx[ip]->flt[i] = 1.;
	}
    }
	
    
    nhelicon_lop(false,false,nd,nd,dfre,rr);
    for (i=0; i < nd; i++) {
	if ( rr[i] > 0.) aa->mis[i] = true;
    }

    for (ip=0; ip < aa->np; ip++) {
	for (i=0; i < aa->hlx[ip]->nh; i++) {
	    aa->hlx[ip]->flt[i] = 0.;
	}
    }
}
