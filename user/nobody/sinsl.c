/* Shaping complex sinusoids */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "sin.h"

static int n,k,n2;
static float amp;
static sf_complex *tmp1, *tmp2;

void sinsl_init(int m         /* data dimensions */, 
		int rect      /* triangle radius */,
		sf_complex z0 /* center frequency */)
/*< initialize >*/
{
    n = m;
    n2 = n+2*rect;

    tmp1 = sf_complexalloc(n2);
    tmp2 = sf_complexalloc(n2);

    k = rect;
    amp = 1.0/rect;

    sinpred_init (z0,n,k);
}

void sinsl_close(void)
/*< free allocated storage >*/
{
    free(tmp1);
    free(tmp2);
}

void sinsl_lop(bool adj, bool add, int nx, int ny, 
	       sf_complex* x, sf_complex* y)
/*< linear operator >*/
{
    int i;
    sf_complex *inp, *out;

    if (nx != n || ny !=n) 
	sf_error("%s: size problem",__FILE__);

    sf_cadjnull(adj,add,nx,ny,x,y);

    if (adj) {
	inp = y;
	out = x;
    } else {
	inp = x;
	out = y;
    }

    for (i=0; i < k; i++) {
	tmp1[i] = sf_cmplx(0.,0.);
	tmp1[i+k+n] = sf_cmplx(0.,0.);
    }
    for (i=0; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H	 
	tmp1[i+k] = amp*inp[i];
#else
	tmp1[i+k] = sf_crmul(inp[i],amp);
#endif
    }

    sinpredicter_lop  (true,  false, n2, n2, tmp2, tmp1);
    sinsubtracter_lop (true,  false, n2, n2, tmp1, tmp2);
    sinsubtracter_lop (false, false, n2, n2, tmp1, tmp2);
    sinpredicter_lop  (false, false, n2, n2, tmp2, tmp1);

    for (i=0; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H	 
	out[i] += amp*tmp1[i+k];
#else
	out[i] = sf_cadd(out[i], sf_crmul(tmp1[i+k],amp));
#endif
    }
}
