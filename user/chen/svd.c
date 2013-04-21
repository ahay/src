/* singular value decomposition */

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


void svd( int n1, int n2, float **a, bool vt,
	float **u, float *s, float **v)
/*< singular value decomposition >*/
{
	int i1, i2;
	int lwork, info;
	float e, *work;

	lwork = n1*n2;
	work = sf_floatalloc(lwork);

	sgesvd_("A", "A", &n1, &n2, a[0], &n1, 
		s, u[0], &n1, v[0], &n2,
		work, &lwork, &info );
	free(work);

	if(vt == false)
	{
		for(i1=0; i1<n2; i1++)
		for(i2=0; i2<i1; i2++)
		{
			e = v[i2][i1];
			v[i2][i1] = v[i1][i2]; 
			v[i1][i2] = e;
		}
	}
}



