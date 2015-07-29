/* frequency response of PWD operator */

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
#include "lpwd.h"
#include "opwd.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, nf, n3, i3;
	sf_complex **buf, c1;
	float d3, o3, radius;
	bool iir, opwd;
	char *interp;

	sf_init(argc, argv);

	out = sf_output("out");

	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getint("order", &nf)) nf=1;
	/* order of PWD */
	if(!sf_getfloat("o3", &o3)) o3=20;
	/* first dip angle */
	if(!sf_getfloat("d3", &d3)) d3=30;
	/* dip angle increment */
	if(!sf_getint("n3", &n3)) n3=1;
	/* number dip angle */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	sf_putint(out, "n1", 2*n1+1);
	sf_putfloat(out, "o1", -0.5);
	sf_putfloat(out, "d1", 0.5/n1);
	sf_putint(out, "n2", 2*n1+1);
	sf_putfloat(out, "o2", -0.5);
	sf_putfloat(out, "d2", 0.5/n1);
	sf_putint(out, "n3", n3);
	sf_putfloat(out, "d3", d3);
	sf_putfloat(out, "o3", o3);

	buf = sf_complexalloc2(2*n1+1, 2*n1+1);
	sf_settype(out,SF_COMPLEX);

	if(!sf_getbool("iir", &iir)) iir=false;
	/* y: iir; n: fir */

	if(!sf_getbool("opwd", &opwd)) opwd=true;
	/* y: circle interpolating pwd; n: line interpolating pwd */
	if(!sf_getfloat("radius", &radius)) radius=1.0;
	/* radius for circle interpolating pwd */

	if(opwd==true)
	{
		opwd_init(nf, nf, interp, 1.0);
		for(i3=0; i3<n3; i3++)
		{
			c1 = sf_cmplx(0, (d3*i3+o3)/180*SF_PI);
#ifdef SF_HAS_COMPLEX_H
			opwd_freq(radius*cexpf(c1), n1, buf, iir);
#else
			opwd_freq(sf_crmul(cexpf(c1),radius), n1, buf, iir);
#endif
			sf_complexwrite(buf[0], (2*n1+1)*(2*n1+1), out);
		}
		opwd_close();
	}else{
		lpwd_init(nf, nf, 0, 0, interp);
		for(i3=0; i3<n3; i3++)
		{
			lpwd_freq((d3*i3+o3)/180*SF_PI, n1, buf, iir);
			sf_complexwrite(buf[0], (2*n1+1)*(2*n1+1), out);
		}
		lpwd_close();
	}
	free(buf[0]);
	free(buf);

	exit(0);
}


