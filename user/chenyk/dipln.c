/* large dip calculation via non-stationary regularization  */
/*
  Copyright (C) 2012 The University of Texas at Austin
  
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
#include "dipl.h"
#include "divnn.h"

static int n1, n2, nf1, nf2;
static float **u1, **u2, **u3, **u4, **u5, r, p0, eps;
static sf_complex **p;
static bool verb, use_divn;

#define divn(a, b)  (a*b/(b*b+10E-15))
/*< initialize >*/

static void odip_verb(int it, int n, float *u)
/*< initialize >*/
{
    int i1;
    double norm;

    if(verb)
    {
	for(i1=0, norm=0.0; i1<n; i1++)
	    norm += (u[i1]*u[i1]);
	sf_warning("res1 %d %g", it, sqrtf(norm/n));
    }
}

static int nn[3];

/*for non-stationary regularization*/
void odipn_init(char* interp, int mf1, int mf2, float rad,
	       int m1, int m2, int *nbox, float **rct, int **sft , int niter, float dip0, float eps1, bool vb)
/*< initialize >*/
{
/* int *nbox: triangle radius [ndim] */
/* float **rct: triangle lengths [ndim][nd] */
/* int **sft: triangle shifts [ndim][nd] */	 

    int n;
    nf1 = mf1;
    nf2 = mf2;
    n1 = m1;
    n2 = m2;
    verb = vb;
    eps=eps1;

    u1 = sf_floatalloc2(n1, n2);
    u2 = sf_floatalloc2(n1, n2);
    u3 = sf_floatalloc2(n1, n2);
    u4 = sf_floatalloc2(n1, n2);
    u5 = sf_floatalloc2(n1, n2);
    p = sf_complexalloc2(n1, n2);

    r=rad;
    p0 = dip0;
	
    opwd_init(nf1, nf2, interp, r);
    if(nbox[0]>0 && nbox[1]>0)
    {
	n = n1*n2;
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = 1;
// 	sf_divn_init (2, n, nn, rect, niter, false);
	divnn_init (3, n, nbox, nn, rct, sft, niter, verb);
// 	sf_divn_init (2, n, nn, nbox, niter, false);
	use_divn=true;
    }else 	use_divn=false;
}

void odipn_close()
/*< release memory >*/
{
    free(u1[0]);
    free(u2[0]);
    free(u3[0]);
    free(u4[0]);
    free(u5[0]);
    free(u1);
    free(u2);
    free(u3);
    free(u4);
    free(u5);
    free(p[0]);
    free(p);
    opwd_close();
    if(use_divn)	divnn_close();
}


static void odipn_verb(int it, int n, float *u)
/*< initialize >*/
{
    int i1;
    double norm;

    if(verb)
    {
	for(i1=0, norm=0.0; i1<n; i1++)
	    norm += (u[i1]*u[i1]);
	sf_warning("res1 %d %g", it, sqrtf(norm/n));
    }
}

void odipn(float **in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
    int it, i1;
    double  c1;
    sf_complex dip0;

#ifdef SF_HAS_COMPLEX_H
    dip0 = r*cexpf(sf_cmplx(0.0, p0));
#else
    dip0 = sf_crmul(cexpf(sf_cmplx(0.0, p0)),r);
#endif

    for(i1=0; i1<n1*n2; i1++)
    {
	dip[0][i1] = p0;
	p[0][i1] = dip0;
    }

    for (it=0; it<nit; it++)
    {
	opwd(n1, n2, in, p, u1);
	opwdpd(n1, n2, in, p, u2, 0);
	opwdpd(n1, n2, in, p, u3, 1);

	odip_verb(it+1, n1*n2, u1[0]);
	for(i1=0, c1=0.0; i1<n1*n2; i1++)
	{
	    u4[0][i1] = u2[0][i1]*crealf(p[0][i1]) - u3[0][i1]*cimagf(p[0][i1]);
	    c1 += (u4[0][i1]*u4[0][i1]);
	}
	c1=sqrtf(c1/(n1*n2));
	for(i1=0; i1<n1*n2; i1++)
	{
	    u1[0][i1] /= c1;
	    u4[0][i1] /= c1;
	}

	if(use_divn)
	{
	    divnne(u1[0], u4[0], u5[0], eps);
	}else{
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u5[0][i1] = divn(u1[0][i1], u4[0][i1]);
	    }
	}
	for(i1=0; i1<n1*n2; i1++)
	{
	    dip[0][i1] -= eta * u5[0][i1];
#ifdef SF_HAS_COMPLEX_H
	    p[0][i1] = r*cexpf(sf_cmplx(0, dip[0][i1]));
#else
	    p[0][i1] = sf_crmul(cexpf(sf_cmplx(0, dip[0][i1])),r);
#endif
	}

    }
}

void oslopen(float **in, float **dip, int nit, float eta)
/*< omnidirectional slope estimation >*/
{
    int it, i1;
    double  s1, c1;
    sf_complex dip0;

#ifdef SF_HAS_COMPLEX_H
    dip0 = r*cexpf(sf_cmplx(0, p0));
#else
    dip0 = sf_crmul(cexpf(sf_cmplx(0, p0)),r);
#endif
    for(i1=0; i1<n1*n2; i1++)
	p[0][i1] = dip0;

    for (it=0; it<nit; it++)
    {
	opwd(n1, n2, in, p, u1);
	opwdpd(n1, n2, in, p, u2, 0);
	opwdpd(n1, n2, in, p, u3, 1);

	odip_verb(it+1, n1*n2, u1[0]);

	if(use_divn)
	{
	    for(i1=0, c1=0.0, s1=0.0; i1<n1*n2; i1++)
	    {
		c1 += (u2[0][i1]*u2[0][i1]);
		s1 += (u3[0][i1]*u3[0][i1]);
	    }
	    c1=sqrtf(c1/(n1*n2));
	    s1=sqrtf(s1/(n1*n2));
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u1[0][i1] /= c1;
		u2[0][i1] /= c1;
	    }
	    divnne (u1[0], u2[0], u4[0], eps);
	    
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u1[0][i1] *= c1/s1;
		u3[0][i1] /= s1;
	    }
	    divnne (u1[0], u3[0], u5[0], eps);
	}else{
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u4[0][i1] = divn(u1[0][i1], u2[0][i1]);
		u5[0][i1] = divn(u1[0][i1], u3[0][i1]);
	    }
	}
	for(i1=0; i1<n1*n2; i1++)
	{
#ifdef SF_HAS_COMPLEX_H
	    p[0][i1] -= eta * sf_cmplx(u5[0][i1], u4[0][i1]);
	    p[0][i1] = p[0][i1]*r/(cabsf(p[0][i1])+ 1E-15);
#else
	    p[0][i1] = sf_cadd(p[0][i1],sf_crmul(sf_cmplx(u5[0][i1], u4[0][i1]),-eta));
	    p[0][i1] = sf_crmul(p[0][i1],r/(cabsf(p[0][i1])+ 1E-15));
#endif
	}

    }
    for(i1=0; i1<n1*n2; i1++)
	dip[0][i1] = atan(tan(cargf(p[0][i1])));
}



