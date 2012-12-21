/* omnidirectional dip estimation  */

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
#include "fbdip.h"


int main(int argc, char*argv[])
{
	sf_file in, out;
	int n1, n2, n3, m1, m2, rect[2], niter, liter;
	int i3;
	bool verb;
	float ****wav, **dip, radius, eta, idip;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &m1)) sf_error("No m1= in input");
	if (!sf_histint(in, "n4", &m2)) sf_error("No m2= in input");
	n3 = sf_leftsize(in, 4);
	sf_unshiftdim2(in,out,3);

	if (!sf_getint("rect1",&rect[0])) rect[0]=0;
	/* dip smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=0;
	/* dip smoothness on 2nd axis */
	
	if (!sf_getint("niter",&niter)) niter=5;
	/* number of iterations */
	if (!sf_getint("liter",&liter)) liter=20;
	/* number of linear iterations */
	if (!sf_getfloat("radius", &radius)) radius = 1.0;
	/* interpolating radius for opwd */
	if (!sf_getfloat("eta", &eta)) eta = 0.75;
	/* steps for iteration */
	if (!sf_getfloat("idip", &idip)) idip = 0.0;
	/* initial dip */
	if (!sf_getbool("verb", &verb)) verb = false;
	/* verbosity flag */

	wav = sf_floatalloc4(n1, n2, m1, m2);
	dip = sf_floatalloc2(n1,n2);

	/* initialize dip estimation */
	fbdip_init(radius, n1, n2, rect, liter, idip, verb);


	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0][0][0], n1*n2*m1*m2, in);
		fbdip(m1, m2, wav, dip, niter, eta);
		sf_floatwrite(dip[0], n1*n2, out);
	}

	fbdip_close();
	free(dip[0]);
	free(dip);
	free(wav[0][0][0]);
	free(wav[0][0]);
	free(wav[0]);
	free(wav);
	return 0;
}



