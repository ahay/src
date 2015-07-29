/* frequency response of linear phase approximators */

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
#include "lphpoly.h"
#include "dsp.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, n2, nf, m, n, i1, i2;
	sf_complex z1, **bf;
	float o1, d1, o2, d2, **c, *b1, f, delay;
	bool iir;
	char* interp;

	sf_init(argc, argv);

	out = sf_output("out");
	sf_settype(out,SF_COMPLEX);

	if(!sf_getint("m", &m)) m=1;
	/* b[-m, ... ,n] */
	if(!sf_getint("n", &n)) n=1;
	/* b[-m, ... ,n] */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	nf = m+n+1;

	if(!sf_getbool("iir", &iir)) iir=true;
	/* y:iir,  n:fir */
	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getfloat("o2", &o2)) o2=0.1;
	/* first phase shift */
	if(!sf_getfloat("d2", &d2)) d2=0.3;
	/* phase shift increment */
	if(!sf_getint("n2", &n2)) n2=1;
	/* number of phase shift */

	o1 = 0.0;
	d1 = 0.5/(n1-1);
	sf_putint(out, "n1", n1);
	sf_putfloat(out, "o1", o1);
	sf_putfloat(out, "d1", d1);
	sf_putint(out, "n2", n2);
	sf_putfloat(out, "o2", o2);
	sf_putfloat(out, "d2", d2);

	bf = sf_complexalloc2(n1, n2);
	b1 = sf_floatalloc(nf);

	c = lphpoly(m, n, interp);
	if(c==NULL) sf_error("interp=%s not correct", interp);
	for(i2=0; i2<n2; i2++)
	{
		delay = o2+d2*i2;
		if(iir) delay /= 2.0;
		lphpoly_coef(nf-1, c, delay, b1, nf-1, false);
		for(i1=0; i1<n1; i1++)
		{
			f = (o1+d1*i1);
			z1 = fir_freq(-m, n, b1+m, f);
#ifdef SF_HAS_COMPLEX_H
			if(iir) z1 = z1*z1/(z1*conjf(z1)+0.000001);
#else
			if(iir) z1 = sf_crmul(sf_cmul(z1,z1),1.0/(crealf(sf_cmul(z1,conjf(z1)))+0.000001));
#endif
			bf[i2][i1] = z1;
		}
	}

	sf_complexwrite(bf[0], n1*n2, out);
	free(bf[0]);
	free(bf);
	free(c[0]);
	free(c);
	free(b1);
	exit(0);
}




