/* slope scan */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

float stack(int n1, int n2, float **p)
{
	int i1, i2;
	float f1, f2, t1, f3;

	f2 = 0.0; f3 = 0.0;
	for(i1=0; i1<n1; i1++)
	{
		f1 = 0.0;
		for(i2=0; i2<n2; i2++)
		{
			t1 = p[i2][i1];
			f1 += t1;
			f2 += t1*t1;
		}
		f3 += f1*f1;
	}
	f1 = f2*n2;
	if(f1 <= 0.0) return 0.0;
	else return  (f3/f1);
}

static	float **p, *pv;
static int r1, r2, nw;

void slopescan_init(int rc1, int rc2, int nk, int stk)
/*< initiliaze >*/
{
	r1 = rc1;
	r2 = rc2;
	nw = nk;

	p = sf_floatalloc2((2*r1+1), (2*r2+1));
	pv = sf_floatalloc(2*nw+1);
}

void slopescan_close()
/*< release memory >*/
{
	free(pv);
	free(p[0]);
	free(p);
}

void slopescan(int n1, int n2, float **d1, float **d2)
/*< slope scan >*/
{
	int i1, i2, k1, j1, j2, j3, j4, mx;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		for(k1=-nw, mx=-nw; k1<nw; k1++)
		{
			for(j1=-r1; j1<r1; j1++)
			for(j2=-r2; j2<r2; j2++)
			{
				j3 = i1+j2*k1+j1;
				j4 = i2+j2;
				if(j3<0 || j4<0 || j3>=n1 || j4>=n2) p[j2+r2][j1+r1] = 0.0;
				else p[j2+r2][j1+r1] = d1[j4][j3];
			}
			pv[k1+nw] = stack(2*r1+1, 2*r2+1, p);
			if(pv[k1+nw] > pv[mx+nw]) mx = k1;
		}
		d2[i2][i1] = mx;
	}
}


