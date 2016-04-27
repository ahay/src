/* Update model with search direction and step length in FWI*/
/*
 Copyright (C) 2014 University of Texas at Austin
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include "optimization.h"

int main(int argc, char* argv[])
{
	int n1, n2;
	float **m0, **m1, **d, alpha0, max;
	sf_file in, out, dir, alpha;

	sf_init(argc, argv);

	in=sf_input("in");
	out=sf_output("out");
	dir=sf_input("direction");
	alpha=sf_input("alpha");

	if(!sf_getfloat("max", &max)) max=1.;
	if(!sf_histint(in, "n1", &n1)) sf_error("No n1 in input.");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2 in input.");

	m0=sf_floatalloc2(n1, n2);
	m1=sf_floatalloc2(n1, n2);
	d=sf_floatalloc2(n1, n2);
	sf_floatread(m0[0], n1*n2, in);
	sf_floatread(d[0], n1*n2, dir);
	sf_floatread(&alpha0, 1, alpha);

	update_model_fwi(m0,m1,d,alpha0,max,n1,n2);

	sf_floatwrite(m1[0], n1*n2, out);

	exit(0);
}
