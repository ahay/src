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

#include "expsignoi3.h"
#include "expont4.h"

static int n1,n2,n3,n123;
static float eps, *a, *b, *c, *d, *tmp, *aq, *bq, *cq, *dq;

void expsignoi3_init (int m1,int m2,int m3 /* data size */, 
		      float eps1    /* signal/noise scaling */, 
		      float **aa    /* frequency filter [8][m1*m2] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    n3 = m3;
    n123 = n1*n2*n3;
    eps = eps1;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    aq = aa[4];
    bq = aa[5];
    cq = aa[6];
    dq = aa[7];

    tmp = sf_floatalloc(n123);
}

void expsignoi3_close(void)
/*< free allocated storage >*/
{
    free (tmp);
}

void expsignoi3_lop (bool adj, bool add, int ns, int nd, 
		     float *sig /* signal */, float *dat /* data */)
/*< linear operator >*/
{
    int i;

    if (4*ns != nd) sf_error("%s: wrong size: 4*%d != %d",__FILE__,ns,nd);

    sf_adjnull(adj,add,ns,nd,sig,dat);

    expont4_init (n1, n2, n3, a, b, aq, bq);
    expont4_lop (adj, true, ns, 2*ns, sig,dat);
    expont4_init (n1, n2, n3, c, d, cq, dq);

    if (adj) {
    for (i=0; i < n123; i++) {
	tmp[i] = dat[i+n123+n123] * eps;
    }
    for (i=n123; i < 2*n123; i++) {
        tmp[i] = dat[i+n123+n123] * eps;
    }
    } else {
    for (i=0; i < n123; i++) {
        tmp[i] = sig[i] * eps;
    }
    }

    if (adj) {
	expont4_lop (true, true, ns, 2*ns, sig, tmp);
    } else {
	expont4_lop (false, true, ns, 2*ns, tmp, dat+n123);
    }
}
