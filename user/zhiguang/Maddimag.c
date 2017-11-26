/* Convert large-size (with n3=) real data to complex (by adding zero imaginary part) */
/*
 Copyright (C) 2015 University of Texas at Austin
 
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
	int i3, n3, ix, nx, iz, nz;
	float **a;
	sf_complex **b;
	sf_file in, out;

	sf_init(argc, argv);

	in=sf_input("in");
	out=sf_output("out");
	sf_settype(out, SF_COMPLEX);

	if(!sf_histint(in, "n1", &nz)) sf_error("No n1= in input");
	if(!sf_histint(in, "n2", &nx)) sf_error("No n2= in input");

	n3=sf_leftsize(in, 2);

	a=sf_floatalloc2(nz, nx);
	b=sf_complexalloc2(nz, nx);

	for (i3=0; i3<n3; i3++){
		sf_floatread(a[0], nz*nx, in);

#pragma omp parallel for private(ix, iz)
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				b[ix][iz]=sf_cmplx(a[ix][iz],0.);
			}
		}

		sf_complexwrite(b[0], nz*nx, out);
	}

	exit(0);
}
