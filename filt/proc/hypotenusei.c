#include <math.h>

#include <rsf.h>

#include "hypotenusei.h"

static int nt, *iz;

void hypotenusei_init(int nt1)
{
    nt = nt1;
    iz = sf_intalloc(nt);
}

void hypotenusei_set(float t0, float dt, float xs)
{
    int it;
    float t, zsquared;

    for (it=0; it < nt; it++) {  
	t = t0 + dt*it;
        zsquared =  t * t - xs * xs;
	iz[it] = ( zsquared >= 0.)? 0.5 + (sqrtf( zsquared) - t0) /dt: -1;
    }
}

void hypotenusei_close(void)
{
    free(iz);
}

void hypotenusei_lop(bool adj, bool add, int n1, int n2, float *zz, float *tt)
{
    int  it;
    
    sf_adjnull(adj,add,n1,n2,zz,tt);

    for (it=0; it < nt; it++) {  
	if (iz[it] < 0) continue;

	if (adj) 
	    zz[iz[it]] +=  tt[it];
	else 
	    tt[it] +=  zz[iz[it]];
    }
}

