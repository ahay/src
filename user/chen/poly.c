/* polynomial operators */

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

/* p(x) = c[0] + c[1]*x + c[2]*x^2 + c[n]*x^n */

#include <rsf.h>

static double prod_n_m(double *r, int n, int m)
{
    if( n<m || n<=0) return 0.0;
    if( m == 0 ) return 1.0;
    if( m == 1 ) 
	return (r[0]+prod_n_m(r+1, n-1, 1));
    return (r[0]*prod_n_m(r+1, n-1, m-1) +
	    prod_n_m(r+1, n-1, m));
}

void poly(int n, float *r)
/*< polynomial root to coefficients 
  [1.0, 2.0, 1.0] => [2.0, 3.0, 1.0]
  >*/
{
    double *p;
    int i1;

    p = sf_alloc(n, sizeof(double));
    for(i1=0; i1<n; i1++)	p[i1] = -r[i1];

    for(i1=0; i1<=n; i1++)
    {
	r[i1]=r[n]*prod_n_m(p, n, n-i1);
    }

    free(p);
}

sf_complex poly_dval(int n, float *c, sf_complex x)
/*< polynomial derivative value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(0.0,0.0);
    for(i=n; i>=1; i--)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*x + c[i]*i;
#else
	val = sf_cadd(sf_cmul(val,x),sf_cmplx(c[i]*i,0.0));
#endif
    }
    return val;
}



sf_complex poly_val(int n, float *c, sf_complex x )
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(0.0,0.0);
    for(i=n; i>=0; i--)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*x + c[i];
#else
	val = sf_cadd(sf_cmul(val,x),sf_cmplx(c[i],0.0));
#endif
    }
    return val;
}

sf_complex plyr_val(int n, float *r, sf_complex x)
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(r[n],0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*(x-r[i]);
#else
	val = sf_cmul(val,sf_cadd(x,sf_cmplx(-r[i],0.0)));
#endif
    }
    return val;
}


sf_complex plyr_dval(int n, float *r, sf_complex x)
/*< polynomial derivative value >*/
{
    int i, k;
    sf_complex v1, v2;

    v1 = sf_cmplx(1.0,0.0);
    k = 0;
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 = x - r[i];
#else
	v2 = sf_cadd(x,sf_cmplx(-r[i],0.0));
#endif
	if(crealf(v2) == 0.0f && cimagf(v2) == 0.0f) k++;
#ifdef SF_HAS_COMPLEX_H
	else	v1 *= v2;
#else
	else    v1 = sf_cmul(v1,v2);
#endif
    }
    if(k >= 2) return sf_cmplx(0.0,0.0);
#ifdef SF_HAS_COMPLEX_H
    if(k == 1) return v1*r[n];
#else
    if(k == 1) return sf_crmul(v1,r[n]);
#endif

    v2 = sf_cmplx(0.0,0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 += v1/(x-r[i]);
#else
	v2 = sf_cadd(v2,sf_cdiv(v1,sf_cadd(x,sf_cmplx(-r[i],0.0))));
#endif
    }

#ifdef SF_HAS_COMPLEX_H
    return v2*r[n];
#else
    return sf_crmul(v2,r[n]);
#endif
}


sf_complex cplyr_val(int n, sf_complex *r, sf_complex x)
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = r[n];
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*(x-r[i]);
#else
	val = sf_cmul(val,sf_cadd(x,sf_crmul(r[i],-1.0)));
#endif
    }
    return val;
}

sf_complex cplyr_dval(int n, sf_complex *r, sf_complex x)
/*< polynomial derivative value >*/
{
    int i, k;
    sf_complex v1, v2;

    v1 = sf_cmplx(1.0,0.0);
    k = 0;
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 = x - r[i];
#else
	v2 = sf_cadd(x,sf_crmul(r[i],-1.0));
#endif
	if(crealf(v2) == 0.0f && cimagf(v2) == 0.0f) k++;
#ifdef SF_HAS_COMPLEX_H
	else	v1 *= v2;
#else
	else	v1 = sf_cmul(v1,v2);
#endif
    }
    if(k >= 2) return sf_cmplx(0.0,0.0);
#ifdef SF_HAS_COMPLEX_H
    if(k == 1) return v1*r[n];
#else
    if(k == 1) return sf_cmul(v1,r[n]);
#endif

    v2 = sf_cmplx(0.0,0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 += v1/(x-r[i]);
#else
	v2 = sf_cadd(v2,sf_cdiv(v1,sf_cadd(x,sf_crmul(r[i],-1.0))));
#endif
    }

#ifdef SF_HAS_COMPLEX_H
    return v2*r[n];
#else
    return sf_cmul(v2,r[n]);
#endif
}



