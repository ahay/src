#include <assert.h>

#include <rsf.h>

#include "igrad2.h" 
/* 2-D gradient with adjoint,  r= grad( p) */

static int n1, n2, n12; 

void igrad2_init (int n1_in, int n2_in)
{
    n1 = n1_in; 
    n2 = n2_in;
    n12 = n1*n2;
}

void igrad2_lop (bool adj, bool add, int np, int nr, float* p, float* r)
{
    int i1,i2,i;

    assert (np == n12);
    assert (nr == 2*n12);

    sf_adjnull (adj,add,np,nr,p,r);

    for (i2=0; i2 < n2-1; i2++) {  
	for (i1=0; i1 < n1-1; i1++) {
	    i = i1+i2*n1;
	    if (adj == true) {
		p[i+1]  += r[i]; 
		p[i+n1] += r[i+n12];
		p[i]    -= (r[i] + r[i+n12]);
	    } else {
		r[i]     += (p[i+1]  - p[i]); 
		r[i+n12] += (p[i+n1] - p[i]);
	    }
	}
    }
}
