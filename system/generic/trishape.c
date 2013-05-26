/* 2-D irregular data interpolation using triangulation and shaping regularization. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

#include "list_struct.h"
#include "delaunay.h"

static bool sym;
static int nd, n1, n2, n12;
static float o1,d1, o2,d2;
static float *d, *m;
static Node q=NULL;
static sf_triangle tr1, tr2;

void trishape_init(bool sym1, int nd1, int n11, int n21, 
		   float o11, float o21, float d11, float d21,
		   int rect1, int rect2, int nw,
		   float **xy)
/*< initialize >*/
{
    sym = sym1;
    nd = nd1;
    n1 = n11;
    n2 = n21;
    n12 = n1*n2;

    o1 = o11; d1=d11;
    o2 = o21; d2=d21;

    tr1 = (rect1 > 1)? sf_triangle_init (rect1,n1): NULL;
    tr2 = (rect2 > 1)? sf_triangle_init (rect2,n2): NULL;

    sf_int2_init (xy, o1, o2, d1, d2, n1, n2, sf_lg_int, nw, nd);

    d = sf_floatalloc(nd);
    if (sym) m = sf_floatalloc(n12);
}

void trishape_close(void)
/*< free allocated storage >*/
{
    if (NULL != tr1) sf_triangle_close(tr1);
    if (NULL != tr2) sf_triangle_close(tr2);

    sf_int2_close();
    free(d);
    if (sym) free(m);
}

void trishape_smooth(float *modl)
/*< shaping operator (triangle smoothing) >*/
{
    int i1, i2;
    
    if (NULL != tr1) {
	for (i2=0; i2 < n2; i2++) {
	    sf_smooth2 (tr1, 0, 1, false, false, modl+i2*n1);
	}
    }
    if (NULL != tr2) {
	for (i1=0; i1 < n1; i1++) {
	    sf_smooth2 (tr2, i1, n1, false, false, modl);
	}
    }   
}    

void trishape_back(const float *data, float *modl)
/*< backward operator (triangulation) >*/
{
    int i, i1, i2;

    if (NULL==q) q = AppendNode (0.,0.,0.,EMPTY);

    if (NULL != data) NodeValues(3, nd, data);

    for (i=0; i < n12; i++) {
	i1 = i%n1;
	i2 = i/n1;

	MoveNode (q,  o1+i1*d1,  o2+i2*d2);
	modl[i] = Interpolate (q);
    }    
}

void trishape_forw(const float *modl, float *data)
/*< forward operator (grid interpolation) >*/
{
    sf_int2_lop (false,false,n12,nd,(float*) modl,data);
}

void trishape(int n12, const float *inp, float *out, void *data)
/*< I + S (BF - I) >*/
{
    int i;

    if (sym) {
	/* (B F - I) S */
	for (i=0; i < n12; i++) {
	    m[i] = inp[i];
	}
	trishape_smooth(m);

	trishape_forw(m,d);
	trishape_back(d,out);

	for (i=0; i < n12; i++) {
	    out[i] -= m[i];
	}
    } else {
	/* B F - I  */

	trishape_forw(inp,d);
	trishape_back(d,out);

	for (i=0; i < n12; i++) {
	    out[i] -= inp[i];
	}
    }

    /* S */
    trishape_smooth(out);

    /* + I */
    for (i=0; i < n12; i++) {
	out[i] += inp[i];
    }
}
