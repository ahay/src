/* stack operators  */

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
#include "sinterp.h"

static sinterp intp=NULL;

void stack_init(sinterp interp)
/*< initialize >*/
{
	intp = interp;
}

typedef float (*stack)(float **u, int n1, int n2, float *c, int win);
/* generic stack interface */
/*^*/

float stack_amp(float ** u, int n1, int n2, float *c, int win)
/*< amplitude stack >*/
{
	double t1;
	int i2;

	t1 = 0.0;
	for(i2=0; i2<n2; i2++)
		t1 += intp(u[i2], c[i2], n1);
	return t1;
}

float stack_l1(float ** u, int n1, int n2, float *c, int win)
/*< stack normalized by L1 norm >*/
{
	double t1, t2, t3;
	int i2;

	t2 = 0.0; t3 = 0.0;
	for(i2=0; i2<n2; i2++)
	{
		t1 = intp(u[i2], c[i2], n1);
		t2 += t1;
		t3 += fabs(t1);
	}
	return (t2/t3);
}

float stack_l2(float ** u, int n1, int n2, float *c, int win)
/*< stack normalized by L2 norm >*/
{
	double t1, t2, t3;
	int i2;

	t2 = 0.0; t3 = 0.0;
	for(i2=0; i2<n2; i2++)
	{
		t1 = intp(u[i2], c[i2], n1);
		t2 += t1;
		t3 += t1*t1;
	}
	return (t2*t2/t3);
}

float stack_mxcor(float ** u, int n1, int n2, float *c, int win)
/*< multi cross correlation  >*/
{
	double t1, t2, t3, t4;
	int i1, i2;

	t4 = 0.0;
	for(i1=-win; i1<=win; i1++)
	{
		t2 = 0.0; t3 = 0.0;
		for(i2=0; i2<n2; i2++)
		{
			t1 = intp(u[i2], c[i2]+i1, n1);
			t2 += t1;
			t3 += t1*t1;
		}
		t4 += (t2*t2-t3);
	}
	return t4/2.0;
}

float stack_ncc(float ** u, int n1, int n2, float *c, int win)
/*< normalized cross correlation  >*/
{
	double t1, t2, t3, t4;
	int i1, i2;

	t3 = 0.0; t4 = 0.0;
	for(i1=-win; i1<=win; i1++)
	{
		t2 = 0.0;
		for(i2=0; i2<n2; i2++)
		{
			t1 = intp(u[i2], c[i2]+i1, n1);
			t2 += t1;
			t3 += t1*t1;
		}
		t4 += t2*t2;
	}
	if(t3 == 0.0) return -1.0/(n2-1);
	else return (t4/t3-1)/(n2-1);
}


float stack_semb(float ** u, int n1, int n2, float *c, int win)
/*< semblance  >*/
{
	double t1, t2, t3, t4;
	int i1, i2;

	t3 = 0.0; t4 = 0.0;
	for(i1=-win; i1<=win; i1++)
	{
		t2 = 0.0;
		for(i2=0; i2<n2; i2++)
		{
			t1 = intp(u[i2], c[i2]+i1, n1);
			t2 += t1;
			t3 += t1*t1;
		}
		t4 += t2*t2;
	}
	return t4/(t3*n2+1E-20);
}


