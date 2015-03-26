/* Pad boundary */
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

int main(int argc, char* argv[])
{
	int ix, iz, i3;
	int nx, nz, n3;
	int top, bottom, left, right;
	int padnx, padnz;
	float dx, dz;
	float x0, z0;
	sf_axis ax, az;

	float **input, **output;
	sf_file Fin, Fout;

	sf_init(argc, argv);
	Fin=sf_input("in");
	Fout=sf_output("out");

	if(!sf_getint("top", &top)) top=0;
	if(!sf_getint("bottom", &bottom)) bottom=0;
	if(!sf_getint("left", &left)) left=0;
	if(!sf_getint("right", &right)) right=0;

	ax=sf_iaxa(Fin, 2); nx=sf_n(ax); dx=sf_d(ax); x0=sf_o(ax);
	az=sf_iaxa(Fin, 1); nz=sf_n(az); dz=sf_d(az); z0=sf_o(az);
	n3=sf_leftsize(Fin, 2);

	padnx=left+nx+right;
	padnz=top+nz+bottom;
	x0=x0-left*dx;
	z0=z0-top*dz;

	sf_setn(ax, padnx); sf_seto(ax, x0);
	sf_setn(az, padnz); sf_seto(az, z0);

	sf_oaxa(Fout, ax, 2);
	sf_oaxa(Fout, az, 1);

	input=sf_floatalloc2(nz, nx);
	output=sf_floatalloc2(padnz, padnx);

	for(i3=0; i3<n3; i3++){
		sf_floatread(input[0], nz*nx, Fin);
		
		for(ix=0; ix<nx; ix++)
			for(iz=0; iz<nz; iz++)
				output[ix+left][iz+top]=input[ix][iz];
		
		for(ix=left; ix<left+nx; ix++)
			for(iz=0; iz<top; iz++)
				output[ix][iz]=output[ix][top];
		
		for(ix=left; ix<left+nx; ix++)
			for(iz=0; iz<bottom; iz++)
				output[ix][iz+top+nz]=output[ix][top+nz-1];
		
		for(iz=0; iz<padnz; iz++)
			for(ix=0; ix<left; ix++)
				output[ix][iz]=output[left][iz];
		
		for(iz=0; iz<padnz; iz++)
			for(ix=0; ix<right; ix++)
				output[ix+left+nx][iz]=output[left+nx-1][iz];
		
		sf_floatwrite(output[0], padnx*padnz, Fout);
	}

	return 0;
}
