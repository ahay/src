/* FIR filter */

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
#include "rfir.h"
#include "dsp.h"
#include "dim3axis.h"



int main(int argc,char*argv[])
{
	sf_file in, out, fir;
	int n[3], nf, of, dim, axis, i2, i3;
	float *u1, *u2, *c;
	void *h;

	sf_init(argc,argv);

	in  = sf_input("in");
	out = sf_output("out");

	fir = sf_input("fir");
	if (!sf_histint(fir, "n1", &nf)) sf_error("No n1= in fir file");
	if (!sf_histint(fir, "o1", &of)) sf_error("No o1= in fir file");

//	sf_warning("of=%d, nf=%d\n", of, nf);
	if (!sf_getint("axis",&axis)) axis=1;
	/* apply fir filter on which dimension */
	dim = dim3axis(in, axis, n);
	if(dim < axis) sf_error("dim < axis");

	if(axis == 1)
	{
		u1  = sf_floatalloc(n[1]);
		u2  = sf_floatalloc(n[1]);
	}else{
		u1  = sf_floatalloc(n[0]);
		u2  = sf_floatalloc(n[0]);
	}

	c  = sf_floatalloc(nf);

	sf_floatread(c, nf, fir);
	for(i3=0; i3<n[2]; i3++)
	if(axis == 1)
	{
		sf_floatread(u1, n[1], in);
		firs(of, nf+of-1, c-of, u1, 1, n[1], u2, 1);
		sf_floatwrite(u2, n[1], out);
	}else{
		h = rfir_init(nf, c, n[0]);
		for(i2=0; i2<n[1]; i2++)
		{
			sf_floatread(u1, n[0], in);
			rfir(h, u1);
			if(i2+of>=0) sf_floatwrite(u1, n[0], out);
		}
		for(i2=of; i2<0; i2++)
			sf_floatwrite(u1, n[0], out);
		rfir_close(h);
	}
	free(u1);
	free(u2);
	free(c);

	exit(0);
}

