/* ADaptive Eage Preserving Filter */

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
#include "epfad.h"

int main(int argc, char*argv[])
{
	sf_file in,  out;
	int i3, n1, n2, n3;
	float **u1;
	int rect[3];

	sf_init(argc, argv);

	in = sf_input("in");
	out = sf_output("out");

	if(!sf_getint("rect1", &rect[0])) rect[0]=0;
	if(!sf_getint("rect2", &rect[1])) rect[1]=0;
	if(!sf_getint("rect3", &rect[2])) rect[2]=0;
	/* smoothing half window on the #-th dimension */

	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	if(!sf_histint(in, "n2", &n2)) {n2=1; n3=1; rect[1]=0; rect[2]=0;}
	if(!sf_histint(in, "n3", &n3)) {n3=1; rect[2]=0;}

	u1 = sf_floatalloc2(n1, n2);
	epfad_init(n1, n2, rect);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		epfad(u1);
		sf_floatwrite(u1[0], n1*n2, out);
	}
	epfad_close();
	free(u1[0]);
	free(u1);
	return 0;

}


