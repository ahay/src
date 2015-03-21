/* 2-D Laplacian operator for complex numbers */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "laplac2.h"

static int n1, n2;

void laplac2_init(int m1, int m2)
/*< initialize with data size >*/
{
    n1 = m1;
    n2 = m2;
}

void laplac2(int   np, 
             int   nr, 
             sf_complex *p, 
             sf_complex *r)
/*< linear operator >*/
{
  int i1,i2,j;

  for (i2=0; i2 < n2; i2++) {
    for (i1=0; i1 < n1; i1++) {
      j = i1+i2*n1;
      r[j] = sf_cmplx(0.,0.);
    }
  }

  for (i2=0; i2 < n2; i2++) {
    for (i1=0; i1 < n1; i1++) {
      j = i1+i2*n1;

      if (i1 > 0) {
#ifdef SF_HAS_COMPLEX_H
        r[j] += p[j] - p[j-1];
#else
        r[j] = sf_cadd(r[j],sf_csub(p[j],p[j-1]));
#endif
      }
      if (i1 < n1-1) {
#ifdef SF_HAS_COMPLEX_H
        r[j] += p[j] - p[j+1];
#else
        r[j] = sf_cadd(r[j],sf_csub(p[j],p[j+1]));
#endif
      }

      if (i2 > 0) {
#ifdef SF_HAS_COMPLEX_H
        r[j] += p[j] - p[j-n1];
#else
        r[j] = sf_cadd(r[j],sf_csub(p[j],p[j-n1]));
#endif
      }
      if (i2 < n2-1) {
#ifdef SF_HAS_COMPLEX_H
        r[j] += p[j] - p[j+n1];
#else
        r[j] = sf_cadd(r[j],sf_csub(p[j],p[j+n1]));
#endif
      }
    }
  }
}
