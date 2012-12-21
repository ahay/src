/* 2D dip filter */

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
#include "dipflt.h"


int main(int argc, char*argv[])
{
	sf_file in, out, dip;
	int n1, n2, n3, m1, m2;
	int i3, nf
;
	char *interp, *filt;
	float **u1, **u2, **u3;
	void *h;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");
	dip = sf_input("dip");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type input");
	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (SF_FLOAT != sf_gettype(dip)) sf_error("Need float type dip");
	if (!sf_histint(dip, "n1", &m1)) sf_error("No n1= in dip");
	if (!sf_histint(dip, "n2", &m2)) sf_error("No n2= in dip");

	if(m1 != n1) sf_error("dip.n1=%d and in.n1=%d",m1,n1);
	if(m2 != n2) sf_error("dip.n2=%d and in.n2=%d",m2,n2);
	n3 = sf_leftsize(in, 2);
	
	if ((interp = sf_getstring("interp")) == NULL) interp = "nearest";
	/* interpolation method: [nearest],linear */
	if ((filt = sf_getstring("filt")) == NULL) filt = "median";
	/* filter type: [median],mean */
	if (!sf_getint("nf", &nf)) nf = 1;
	/* filter length */

	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, n2);
	u3 = sf_floatalloc2(n1, n2);

	/* initialize dip estimation */
	h = dipflt_init(nf, n1, n2, filt, interp);
	sf_floatread(u3[0], n1*n2, dip);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		dipflt(h, u3, u1, u2);
		sf_floatwrite(u2[0], n1*n2, out);
	}

	dipflt_close(h);
	free(u1[0]);
	free(u2[0]);
	free(u3[0]);
	free(u1);
	free(u2);
	free(u3);
	return 0;
}



