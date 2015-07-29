

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
#include "opwd.h"

#define n1	3
#define n2	3

int main(int argc, char*argv[])
{
	int i1, i2;
	float **in, **d1, **d2, **u1;
	sf_complex **p1;
	opwd_init(1, 1, "maxflat", 1.0);

/* **c@9
//  0.1666   0.25    0.8333 
//  0.6666   0.00   -0.1666
//  0.1666  -0.25    0.8333 */


	in = sf_floatalloc2(n1, n2);
	u1 = sf_floatalloc2(n1, n2);
	p1 = sf_complexalloc2(n1, n2);
	d1 = sf_floatalloc2(n1, n2);
	d2 = sf_floatalloc2(n1, n2);
/*	fb = sf_floatalloc4(3, 3, n1, n2); */

/*  0.1   0.2   0.3
//  1.1   1.2   1.3
//  2.1   2.2   2.3 */
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		in[i2][i1] = 0.1*(i2+1)+i1;
		p1[i2][i1] = sf_cmplx(1.0,0.0);
	}


	opwd(3,3, in, p1, u1);
/*	u1[1][1] = 0.1; */
	opwdpd(3,3, in, p1, d1, 0);
/*	d1[1][1] = 1.0 */
	opwdpd(3,3, in, p1, d2, 1);
/*	d2[1][1] = 0.1 */


/* large angle model
//  1	2	1.5
//  1	2	1.5
//  1	2	1.5 */
	
	for(i1=0; i1<n1; i1++)
	{
		in[0][i1] = 1.0;
		in[1][i1] = 2.0;
		in[2][i1] = 1.5;
	}
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
		p1[i2][i1] = sf_cmplx(1.0,0.0);


	opwd(3,3, in, p1, u1);
/*	u1[1][1] = 9.0; */
	opwdpd(3,3, in, p1, d1, 0);
/*	d1[1][1] = 0.0 */
	opwdpd(3,3, in, p1, d2, 1);
/*	d2[1][1] = 9.0 */


	return 0;
}


