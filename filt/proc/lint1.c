#include <rsf.h>

#include "lint1.h"

static float o1, d1, *coord;

void lint1_init (float o1_in, float d1_in, float *coord_in)
{
    o1 = o1_in;
    d1 = d1_in;
    coord = coord_in;
}

void lint1_lop  (bool adj, bool add, int nm, int nd, float *mm, float *dd)
{
    int im, id;
    float f, fx, gx;

    sf_adjnull(adj,add,nm,nd,mm,dd);

    for (id=0; id < nd; id++) {
        f = (coord[id]-o1)/d1;     
	im=f;
        if (0 <= im && im < nm-1) {   
	    fx=f-im;   
	    gx=1.-fx;

	    if(adj) {
		mm[im]   +=  gx * dd[id];
		mm[im+1] +=  fx * dd[id];
	    } else {
		dd[id]   +=  gx * mm[im]  +  fx * mm[im+1];
	    }
	}
    }
}

