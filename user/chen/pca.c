/* PCA analysis */

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

int pca_rank(int n0, float *s, float eta)
/*< rank of signal subspace >*/
{
	int i1;
	double e1, e2;

	for(i1=0, e1=0.0; i1<n0; i1++)
		e1 += s[i1];
	e1 *= eta;
	for(i1=0, e2=0.0; i1<n0; i1++)
	{
		if(e2 > e1) return (i1+1);
		e2 += s[i1];
	}
	return i1;
}

void pca_klt(int n1, int n2, int nc, float *s, float **vt, float **a)
/*< KL transform: a[nc, n2] == s[nc, n2] * vt[n2, n2]  >*/
{
	int i1, i2;

	for(i1=0;i1<nc;i1++)
	for(i2=0;i2<n2;i2++)
	{
		a[i2][i1] = s[i1]*vt[i2][i1];
	}
}

void pca_iklt(int n1, int n2, int nc, float **u, float **a, int lda, float **b)
/*< Inverse KL transform: b[n1, n2] == u[nc, n1]^T * a[nc,n2]  >*/
{
	float alpha, beta;

	alpha = 1.0; beta = 0.0;

	sgemm_("N", "N", &n1, &n2, &nc,
		&alpha, u[0], &n1, a[0], &lda, 
		&beta, b[0], &n1);
}


