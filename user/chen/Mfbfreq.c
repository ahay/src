/* frequency response of linear phase filter bank */

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
	int n1, nd, m, n, nf, i1, j1, id, nfd, n1d, kf, k1;
	sf_complex z1, **bf, *buf;
	float o1, d1, **c, f;
	char* interp, tmp[8];

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

	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain is 2*n1+1 */
	if(!sf_getint("nd", &nd)) nd=1;
	/* nd dimensional filter bank, nd should not be large, or you will need to buy a new disk array */

	o1 = -0.5;
	d1 = 0.5/n1;
	for(j1=1; j1<=nd; j1++)
	{
		sprintf(tmp, "n%d", j1%100u);
		sf_putint(out, tmp, 2*n1+1);
		sprintf(tmp, "o%d", j1%100u);
		sf_putfloat(out, tmp, o1);
		sprintf(tmp, "d%d", j1%100u);
		sf_putfloat(out, tmp, d1);
		sprintf(tmp, "n%d", j1+nd);
		sf_putint(out, tmp, nf);
		sprintf(tmp, "o%d", j1+nd);
		sf_putfloat(out, tmp, 0);
		sprintf(tmp, "d%d", j1+nd);
		sf_putfloat(out, tmp, 1);
	}

	n1 = 2*n1+1;
	bf = sf_complexalloc2(n1, nf);

	c = lphpoly(m, n, interp);
	if(c==NULL) sf_error("interp=%s not correct", interp);
	for(j1=0; j1<nf; j1++)
	{
		for(i1=0; i1<n1; i1++)
		{
			f = (o1+d1*i1);
			z1 = fir_freq(-m, n, c[j1]+m, f);
			bf[j1][i1] = z1;
		}
	}

	nfd=1; n1d=1;
	for(i1=0; i1<nd; i1++) { nfd=nf*nfd; n1d=n1*n1d; }

	buf = sf_complexalloc(n1d);
	for(j1=0; j1<nfd; j1++)
	{
		for(i1=0; i1<n1d; i1++)
		{
			z1 = 1; kf=j1; k1=i1;
			for(id=0; id<nd; id++)
			{
				z1 *= bf[kf%nf][k1%n1];
				kf /= nf;
				k1 /= n1;
			}
			buf[i1] = z1;
		}
		sf_complexwrite(buf, n1d, out);
	}

	free(bf[0]);
	free(bf);
	free(buf);
	free(c[0]);
	free(c);
	return 0;
}

