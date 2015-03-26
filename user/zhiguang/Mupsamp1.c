/* 1-D linear interpolation */
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
#include <math.h>

int main(int argc, char* argv[])
{
	int i2, i, j, k;
	int nm, nd, n2, scale;
	float d1, right;
	float *a, *b;
	sf_file in, out;

	sf_init(argc, argv);

	in=sf_input("in");
	out=sf_output("out");

	if(!sf_getint("scale", &scale)) sf_error("Need scale= !");

	if(!sf_histint(in, "n1", &nm)) sf_error("No n1= in input!");
	if(!sf_histfloat(in, "d1", &d1)) sf_error("No d1= in input!");
	
	nd=(nm-1)*scale+1;
	sf_putint(out, "n1", nd);
	sf_putfloat(out, "d1", d1/scale);
	n2=sf_leftsize(in, 1);

	a=sf_floatalloc(nm);
	b=sf_floatalloc(nd);

	for(i2=0; i2<n2; i2++){
		sf_floatread(a, nm, in);

		for(i=0; i<nd-1; i++){
			j=i/scale;
			k=i%scale;

			right=1.0*k/scale;

			b[i]=a[j]*(1.-right)+a[j+1]*right;
		}
		b[nd-1]=a[nm-1];

		sf_floatwrite(b, nd, out);
	}

	free(a); free(b);
	sf_fileclose(in);
	sf_fileclose(out);

	exit(0);
}
