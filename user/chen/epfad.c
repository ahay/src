/* Adaptive Edge-Preserving Filter */

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
#include "recursion.h"

static int n1, n2, rc1, rc2, rc3;
static void *h;
static float *c1, *c2, *c3;

void epfad_app(float *out, float **in, int n12, int nf, int *par)
{
	int i1, i2, j1, j2;
	float vf, vb, ef, eb, *p, **p2;

	if(rc1>0)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		vf = 0.0;
		vb = 0.0;
		p = in[0]+i2*n1+i1;
		for(j1=0; j1<=rc1; j1++)
		{
			if(i1-j1>=0) vf += p[-j1]*c1[j1];
			if(i1+j1<n1) vb += p[j1]*c1[j1];
		}
//		vf /= (rc1+1);	vb /= (rc1+1);
		ef = vf - p[0];	eb = vb - p[0];
		ef = ef*ef; eb = eb*eb;
		p[0] = (eb*vf + ef*vb);
		p[0] /=(eb+ef);
	}

	if(rc2>0)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		vf = 0.0;
		vb = 0.0;
		p = in[0]+i2*n1+i1;
		for(j2=0; j2<=rc2; j2++)
		{
			if(i2-j2>=0) vf += p[-j2*n1];
			if(i2+j2<n1) vb += p[j2*n1];
		}
		vf /= rc2;	vb /= rc2;
		ef = vf - p[0];	eb = vb - p[0];
		ef *= ef; eb *= eb;
		p[0] = (eb*vf + ef*vb);
		p[0] /= (eb+ef);
	}

	if(rc3>0)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		vf = 0.0;
		vb = 0.0;
		p2 = in+rc3;
		for(j2=0; j2<=rc3; j2++)
		{
			vf += p2[j2][i2*n1+i1];
			vb += p2[-j2][i2*n1+i1];
		}
		vf /= rc3;	vb /= rc3;
		ef = vf - p2[0][i2*n1+i1];
		eb = vb - p2[0][i2*n1+i1];
		ef *= ef; eb *= eb;
		out[i2*n1+i1] = (eb*vf + ef*vb);
		out[i2*n1+i1] /=(eb+ef);
	}else
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
		out[i2*n1+i1] = in[0][i2*n1+i1];
}

void epfad_init(int m1, int m2, int *rect)
/*< initialize >*/
{
	int i;
	float aa;

	rc1 = rect[0];
	rc2 = rect[1];
	rc3 = rect[2];
	n1 = m1;
	n2 = m2;

	h = recursion_init(n1*n2, 2*rc3+1, epfad_app, NULL);
	c1 = sf_floatalloc(rc1+1);
	c2 = sf_floatalloc(rc2+1);
	c3 = sf_floatalloc(rc3+1);

	aa = 0.0;
	if(rc1==0) c1[0] = 1.0;
	else{
		for(i=0; i<=rc1; aa+=c1[i], i++)
			c1[i] = 0.5*(1-cos(2*SF_PI*i/rc1));
		for(i=0; i<=rc1; i++) c1[i] /= aa;
	}
	aa = 0.0;
	if(rc2==0) c2[0] = 1.0;
	else{
		for(i=0; i<=rc2; aa+=c2[i], i++)
			c2[i] = 0.5*(1-cos(2*SF_PI*i/rc2));
		for(i=0; i<=rc2; i++) c2[i] /= aa;
	}
	aa = 0.0;
	if(rc3==0) c3[0] = 1.0;
	else{
		for(i=0; i<=rc3; aa+=c3[i], i++)
			c3[i] = 0.5*(1-cos(2*SF_PI*i/rc3));
		for(i=0; i<=rc3; i++) c3[i] /= aa;
	}
}

void epfad_close()
/*< release memory >*/
{
	recursion_close(h);
}

void epfad(float **u1)
/*< adaptive epf >*/
{
	recursion(h, *u1);
}


