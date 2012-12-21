/* dip estimation by linear phase filter bank */

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


static int n1, n2;
static float **u1, **u2, **u3, **u4, **u5, r;
static sf_complex **p, p0;
static bool verb, use_divn;

void fbdip_init(float rad, int m1, int m2, 
	int *rect, int niter, float dip0, bool vb)
/*< initialize >*/
{
	int n, nn[2];
	n1 = m1;
	n2 = m2;
	verb = vb;

	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, n2);
	u3 = sf_floatalloc2(n1, n2);
	u4 = sf_floatalloc2(n1, n2);
	u5 = sf_floatalloc2(n1, n2);
	p = sf_complexalloc2(n1, n2);

	r=rad;
	p0 = rad*cexpf(sf_cmplx(0, dip0));
	
	if(rect[0]>0 && rect[1]>0)
	{
		n = n1*n2;
		nn[0] = n1;
		nn[1] = n2;
		sf_divn_init (2, n, nn, rect, niter, false);
		use_divn=true;
	}else 	use_divn=false;
}

void fbdip_close()
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
	if(use_divn)	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+10E-15))

#define POW(a, b) (b<0?0.0:(b==0?1.0:pow(a,b)))

void fbdip(int n3, int n4, float ****in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
	int it, i1, j1, j2;
	double  norm, s1, c1;

	for(i1=0; i1<n1*n2; i1++)
		p[0][i1] = p0;

	for (it=0; it<nit; it++)
	{
		for(i1=0; i1<n1*n2; i1++)
		{
			s1 = cimagf(p[0][i1]);
			c1 = crealf(p[0][i1]);
			u1[0][i1] = 0.0;
			u2[0][i1] = 0.0;
			u3[0][i1] = 0.0;
			for(j2=0; j2<n4; j2++)
			for(j1=0; j1<n3; j1++)
			{
				if((j1+j2)%2==0) continue;
				u1[0][i1] += 2.0*in[j2][j1][0][i1]*POW(s1, j1)*POW(c1,j2);
				u2[0][i1] += 2.0*in[j2][j1][0][i1]*j1*POW(s1, j1-1)*POW(c1,j2);
				u3[0][i1] += 2.0*in[j2][j1][0][i1]*POW(s1, j1)*j2*POW(c1,j2-1);
			}
		}

		if(verb)
		{
			for(i1=0, norm=0.0; i1<n1*n2; i1++)
				norm += (u1[0][i1]*u1[0][i1]);
			sf_warning("res1 %d %g", it+1, sqrtf(norm/n1/n2));
		}

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
			sf_divn(u1[0], u2[0], u4[0]);
			for(i1=0; i1<n1*n2; i1++)
			{
				u1[0][i1] *= c1/s1;
				u3[0][i1] /= s1;
			}
			sf_divn(u1[0], u3[0], u5[0]);
		}else{
			for(i1=0; i1<n1*n2; i1++)
			{
				u4[0][i1] = divn(u1[0][i1], u2[0][i1]);
				u5[0][i1] = divn(u1[0][i1], u3[0][i1]);
			}
		}
		for(i1=0; i1<n1*n2; i1++)
		{
			p[0][i1] -= eta * sf_cmplx(u5[0][i1], u4[0][i1]);
			p[0][i1] = p[0][i1]*r/(cabsf(p[0][i1])+ 1E-15);
		}

	}
	for(i1=0; i1<n1*n2; i1++)
		dip[0][i1] = atan(tan(cargf(p[0][i1])));
}



