/* slope estimation by stack scan */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include "slopescan.h"

int main(int argc, char*argv[])
{
	sf_file in, out;
	int i3, n1, n2, n3, nw, rc1, rc2;
	float **u1, **u2;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &n3)) n3=1;

	if (!sf_getint("nw",&nw)) nw=5;
	/* window size is 2*nw+1 */
	if (!sf_getint("rect1",&rc1)) rc1=2;
	/* window size on the 1st dimension */
	if (!sf_getint("rect2",&rc2)) rc2=2;
	/* window size on the 2nd dimension */

	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, n2);

	slopescan_init(rc1, rc2, nw, 0);

	for(i3=0; i3<n3; i3++)
	{	
		sf_floatread(u1[0], n1*n2, in);	
		slopescan(n1, n2, u1, u2);
		sf_floatwrite(u2[0], n1*n2, out);
	}
	slopescan_close();

	free(u1[0]);
	free(u2[0]);
	free(u1);
	free(u2);

	return 0;
}



