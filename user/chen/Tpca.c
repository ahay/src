

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include "svd.c"
#include "pca.h"

#define n1	5
#define n2	3
#define n0	(n1<n2? n1:n2)

int main(int argc, char*argv[])
{
	float **pa, a[]={
  -0.739290  ,-0.696132,  -2.194091  , 0.533211,   0.854345,
  -0.189668 ,  0.628536 ,  1.149399 ,  0.076899 , -2.448373,
   1.186423,  -0.422504  ,-0.987137,  -0.318169  ,-0.756615,};

	float **u, *s, **v, **a1, **a2;

	pa = (float**) sf_alloc(n2, sizeof(float*));
	pa[0] =a; pa[1] = a+n1; pa[2] = a+2*n1;

	a1 = sf_floatalloc2(n1, n2);
	a2 = sf_floatalloc2(n1, n2);
	u = sf_floatalloc2(n1, n1);
	v = sf_floatalloc2(n2, n2);
	s = sf_floatalloc(n0);

	svd(n1, n2, pa, true, u, s, v);
/* u =
//  -0.089041  -0.367808  -0.835785   0.361330   0.166410
//  -0.271957   0.197449   0.082320  -0.105002   0.932337
//  -0.677936   0.625607  -0.196808   0.150797  -0.295880
//   0.083157   0.032841   0.388519   0.913068   0.085829
//   0.672011   0.658231  -0.324046   0.044523   0.090248

// s = 3.4945 1.9044 1.4170

// v =
//   0.675654  -0.345676   0.651153
//  -0.736074  -0.365538   0.569717
//   0.041084  -0.864228  -0.501420 */

	pca_klt(n1, n2, n0, s, v, a1);
/* a1 =
//   2.36107  -2.57220   0.14357
//  -0.65831  -0.69613  -1.64584
//   0.92271   0.80731  -0.71054
//   0.00000   0.00000   0.00000
//   0.00000   0.00000   0.00000 */
	pca_iklt(n1, n2, n0, u, a1, n1, a2);
/* a2 = 
//  -0.739290  -0.189668   1.186423
//  -0.696132   0.628536  -0.422504
//  -2.194091   1.149399  -0.987137
//   0.533211   0.076899  -0.318169
//   0.854345  -2.448373  -0.756615 */
	return 0;
}


