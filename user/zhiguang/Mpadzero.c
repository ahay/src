/* Interpolation from a coarser grid to finer grid with zero padded */
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
	int i3, n3, ix, nx, iz, nz;
	int nx2, nz2, scale;
	float dx, dz;
	float **a, **b;
	sf_file in, out;

	sf_init(argc, argv);

	in=sf_input("in");
	out=sf_output("out");

	if(!sf_getint("scale", &scale)) sf_error("Need scale=");

	if(!sf_histint(in, "n1", &nz)) sf_error("No n1= in input");
	if(!sf_histfloat(in, "d1", &dz)) sf_error("No d1= in input");
	if(!sf_histint(in, "n2", &nx)) sf_error("No n2= in input");
	if(!sf_histfloat(in, "d2", &dx)) sf_error("No d2= in input");

	nx2=(nx-1)*scale+1;
	nz2=(nz-1)*scale+1;
	n3=sf_leftsize(in, 2);

	sf_putint(out, "n1", nz2);
	sf_putfloat(out, "d1", dz/scale);
	sf_putint(out, "n2", nx2);
	sf_putfloat(out, "d2", dx/scale);

	a=sf_floatalloc2(nz, nx);
	b=sf_floatalloc2(nz2, nx2);

	for (i3=0; i3<n3; i3++){
		sf_floatread(a[0], nz*nx, in);

#pragma omp parallel for private(ix, iz)
		for(ix=0; ix<nx2; ix++){
			for(iz=0; iz<nz2; iz++){
				b[ix][iz]=0.;
			}
		}

#pragma omp parallel for private(ix, iz)
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				b[ix*scale][iz*scale]=a[ix][iz];
			}
		}

		sf_floatwrite(b[0], nz2*nx2, out);
	}

	exit(0);
}
