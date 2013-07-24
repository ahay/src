/* 1D polynomial coefficient maxflat filter */

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

typedef struct tag_pcmf
{
	int n;
	double **c;	
}pcmf;
/*    n   2n
//   SUM  SUM  c_ij * Z^ip^j
//   i=-n j=0 */

static double prod_n_m(double f0,double df,int n, int m)
{
	if( n<m || n<=0) return 0.0;
	if(m==1) return (f0*n+df/2*n*(n-1));
	return (f0*prod_n_m(f0+df, df, n-1, m-1) +
			prod_n_m(f0+df, df, n-1, m));
}

void* pcmf_init(int n)
/*< initialize >*/
{
	pcmf *p;
	int m1, m2, k, i, j, m3;
	double d1, d2, d3, d4;

	p = (pcmf*) sf_alloc(1,sizeof(pcmf));
	p->n = n;
	k = 2*n+1;
	p->c = (double**)sf_alloc(k, sizeof(double*));
	p->c[0] = (double*)sf_alloc(k*k, sizeof(double));
	for(i=1; i<k; i++) p->c[i] = p->c[i-1] + k;

	for(k=-n; k<=n; k++)
	{
		m1 = n-k;
		m2 = n+k;
		d1 = 1.0;
		for(i=m1;i>1;i--) d1 /= i;
		for(i=m2;i>1;i--) d1 /= i;
		for(i=0; i<=2*n; i++)
		{
			d2 = 0.0;
			m3=2*n-i;
			for(j=0;j<=m3;j++)
			{
				if((m3-j <= m1) && (j <= m2))
				{
					d3 = prod_n_m(-2*n, 1.0, m1, m3-j);
					d4 = prod_n_m(-2*n, 1.0, m2, j);
					if(d3==0.0) d3=1.0;
					if(d4==0.0) d4=1.0;
					if((m2-j) % 2 == 0) d2 += d3*d4;
					else 	d2 -= d3*d4;
				}
			}
			p->c[k+n][i]=d1*d2;
		}
	}
	return p;
}

void pcmf_filt_1d(void *h, double ag, double *b)
/*< return the 1d filter of angle ag >*/
{
	pcmf *pp;
	int k,i;
	double val, pn, p;

	pp = (pcmf*)h;
	p = tan(ag);
	for(k=0; k<=2*pp->n; k++)
	{
		pn  = 1.0;
		val = 0.0;
		for(i=0; i<=2*pp->n; i++, pn*=p)
		{
			val += pp->c[k][i]*pn;
		}
		b[k] = val;
	}
}

void pcmf_der_1d(void *h, double ag, double *b)
/*< derivatives w.r.t p of the filter >*/
{
	pcmf *pp;
	int k,i;
	double val, pn, p;

	pp = (pcmf*)h;
	p  = tan(ag);

	for(k=0; k<=2*pp->n; k++)
	{
		pn  = 1.0;
		val = 0.0;
		for(i=1; i<=2*pp->n; i++, pn*=p)
		{
			val += pp->c[k][i]*pn*i;
		}
		b[k] = val*(p*p+1);
	}
}

void pcmf_filt_2d(void *h, double ag, double **b)
/*< return the 2d filter of p >*/
{
	pcmf *pp;
	int k1, k2, i;
	double b1, b2, pn1, pn2, p1, p2;

	pp = (pcmf*)h;
	p1 = sin(ag);
	p2 = cos(ag);

	for(k1=0; k1<=2*pp->n; k1++)
	{
		for(k2=0; k2<=2*pp->n; k2++)
		{
			pn1 = 1.0;
			pn2 = 1.0;
			b1  = 0.0;
			b2  = 0.0;
			for(i=0; i<=2*pp->n; i++, pn1*=p1, pn2*=p2)
			{
				b1 += pp->c[k1][i]*pn1;
				b2 += pp->c[k2][i]*pn2;
			}
			b[k1][k2] = b1*b2;
		}
	}
}


void pcmf_der_2d(void *h, double ag, double **b)
/*< derivatives w.r.t p of the filter >*/
{
	pcmf *pp;
	int k1, k2, i;
	double b1, b2, d1, d2,  pn1, pn2, p1, p2;

	pp = (pcmf*)h;
	p1 = sin(ag);
	p2 = cos(ag);

	for(k1=0; k1<=2*pp->n; k1++)
	{
		for(k2=0; k2<=2*pp->n; k2++)
		{
			pn1 = 1.0;	d1  = 0.0;
			pn2 = 1.0;	d2  = 0.0;
			b1 = pp->c[k1][0];
			b2 = pp->c[k2][0];
			for(i=1; i<=2*pp->n; i++, pn1*=p1, pn2*=p2)
			{
				b1 += pp->c[k1][i]*pn1*p1;
				b2 += pp->c[k2][i]*pn2*p2;
				d1 += pp->c[k1][i]*pn1*i;
				d2 += pp->c[k2][i]*pn2*i;
			}
			b[k1][k2] = -b1*d2*p1+b2*d1*p2;
		}
	}
}




void pcmf_close(void *h)
/*< release the memory >*/
{
	free(((pcmf*)h)->c[0]);
	free(((pcmf*)h)->c);
	free(h);
}

