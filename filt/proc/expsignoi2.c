/* Frequency-based signal-noise separation, 2-D */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/

#include "expsignoi2.h"
#include "expont.h"

static int n1,n2,n12;
static float eps, *a, *b, *c, *d, *tmp;

void expsignoi2_init (int m1,int m2 /* data size */, 
		      float eps1    /* signal/noise scaling */, 
		      float **aa    /* frequency filter [4][m1*m2] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    n12 = n1*n2;
    eps = eps1;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    tmp = sf_floatalloc(n12);
}

void expsignoi2_close(void)
/*< free allocated storage >*/
{
    free (tmp);
}

void expsignoi2_lop (bool adj, bool add, int ns, int nd, 
		     float *sig /* signal */, float *dat /* data */)
/*< linear operator >*/
{
    int i;

    if (2*ns != nd) sf_error("%s: wrong size: 2*%d != %d",__FILE__,ns,nd);

    sf_adjnull(adj,add,ns,nd,sig,dat);

    expont_init (n1, n2, a, b);
    expont_lop (adj, true, ns, ns, sig, dat);
    expont_init (n1, n2, c, d);

    for (i=0; i < n12; i++) {
	tmp[i] = adj? dat[i+n12] * eps: sig[i] * eps;
    }

    if (adj) {
	expont_lop (true, true, ns, ns, sig, tmp);
    } else {
	expont_lop (false, true, ns, ns, tmp, dat+n12);
    }
}
