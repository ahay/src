/* Bilateral filter */

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
#include "bilateral.h"


int main(int argc, char*argv[])
{
	sf_file in,  out;
	int i3, n1, n2, n3;
	float **u1, sigmar, sigmad;
	bool twod;
	void *h;

	sf_init(argc, argv);

	in = sf_input("in");
	out = sf_output("out");

	if(!sf_getbool("twod", &twod)) twod=false;
	/* y, 2D smoothing */
	if (!sf_getfloat("sigmar", &sigmar)) sigmar=1.0;
	/* sigma for range filter */
	if (!sf_getfloat("sigmad", &sigmad)) sigmad=1.0;
	/* sigma for domain filter */


	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	if(!sf_histint(in, "n2", &n2)) {n2=1; twod=false;}
	if(twod == false) n3 = sf_leftsize(in, 1);
	else n3 = sf_leftsize(in, 2);

	u1 = sf_floatalloc2(n1, n2);
	h = bilateral_init(n1, n2);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		if(sigmad > 0) bilateral_df(h, sigmad);
		if(sigmar > 0) bilateral_rf(h, sigmar, u1[0]);
		bilateral_norm(h);
		bilateral(h, u1[0]);
		sf_floatwrite(u1[0], n1*n2, out);
	}

	free(u1[0]);
	free(u1);
	return 0;

}


