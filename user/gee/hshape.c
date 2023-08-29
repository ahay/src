/* Helical shaping */
/*
  Copyright (C) 2010 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

#include "hshape.h"

static int nx, nr;
static float eps, *tmp1, *tmp2;
static sf_filter pef;

void hshape_init(int ndat      /* data size */,
		 int rect      /* shaping radius */,
		 sf_filter flt /* helical filter */,
		 float eps1    /* regularization */)
/*< initialize >*/
{
    nx = ndat;
    nr = rect;
    eps = eps1*eps1;
    pef = flt;

    sf_helicon_init (pef);

    tmp1 = sf_floatalloc(nx);
    tmp2 = sf_floatalloc(nx);
}

void hshape_close(void)
/*< free allocated storage >*/
{
    free(tmp1);
    free(tmp2);
}

void hshape_lop(bool adj, bool add, int ninp, int nout, float *inp, float *out)
/*< linear operator >*/
{
    int i;
    
    if (nx != ninp || nx != nout) sf_error("wrong dimensions");
    
    sf_adjnull(adj,add,ninp,nout,inp,out);

    for (i=0; i < nx; i++) {
	if (adj) {
	    tmp1[i] = out[i];
	} else {
	    tmp1[i] = inp[i];
	}
    }
    
    hshape(tmp1);

    for (i=0; i < nx; i++) {
	if (adj) {
	    inp[i] += tmp1[i];
	} else {
	    out[i] += tmp1[i];
	}
    }
}

void hshape(float *data)
/*< shaping in place >*/
{
    int i, ir;
    float r;

    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = -0.5*eps/(r*r);

	for (i=0; i < nx; i++) {
	    tmp1[i] = r*data[i];
	}
	
	sf_helicon_lop(false,false,nx,nx,tmp1,tmp2);
	sf_helicon_lop(true,true,nx,nx,data,tmp2);
    }
}

    
    
