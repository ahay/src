#include <rsf.h>

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


typedef struct tag_fd4
{
	int n1;
	int n2;
	double c0;
	double c11;
	double c21;
	double c12;
	double c22;
}fd4;

void * sf_fd4_init(int n1, int n2, float d1, float d2)
/*< initialize >*/
{
	double t;
	fd4 *p;

	p = (fd4*) malloc(sizeof(fd4));
	if(p==NULL) return NULL;

	t = d1;
	t = 1.0/(t*t);
	p->c11 = 4.0*t/3.0;
	p->c12=  -t/12.0;

	t = d2;
	t = 1.0/(t*t);
	p->c21 = 4.0*t/3.0;
	p->c22=  -t/12.0;

	p->c0  = -2.0 * (p->c11 + p->c12 + p->c21 + p->c22);
	p->n1=n1;
	p->n2=n2;
	return (void *)p;
}

void sf_fd4_laplacian(void * h, float **uin, float **uout)
/*< Laplacian operator, 4th-order finite-difference >*/
{
	int i1, i2;
	fd4 *p;
	p = (fd4*)h;
	for (i2=2; i2 < p->n2-2; i2++) 
	{
		for (i1=2; i1 < p->n1-2; i1++) 
		{
			uout[i2][i1] = 
			p->c11*(uin[i2][i1-1]+uin[i2][i1+1]) +
			p->c12*(uin[i2][i1-2]+uin[i2][i1+2]) +
			p->c21*(uin[i2-1][i1]+uin[i2+1][i1]) +
			p->c22*(uin[i2-2][i1]+uin[i2+2][i1]) +
			p->c0*uin[i2][i1];
		}
	}
}

void sf_fd4_close(void * h)
/*< release memory >*/
{
	free(h);
}
