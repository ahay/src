/* Taper with a trapezoidal filter along the first axis */
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

int main(int argc, char *argv[])
{
	int ix, it, i3, nx, nt, n3, range;
	float **dd, *taper;
	sf_file Fdat, Fnew;

	sf_init(argc, argv);

	Fdat=sf_input("in");
	Fnew=sf_output("out");
	if(!sf_getint("length", &range)) range=15; 
	/* tapering length */

	if(!sf_histint(Fdat, "n1", &nt)) sf_error("No n1= in input");
	if(!sf_histint(Fdat, "n2", &nx)) sf_error("No n2= in input");
	n3=sf_leftsize(Fdat, 2);

	dd=sf_floatalloc2(nt, nx);
	taper=sf_floatalloc(nx);
	for(ix=0; ix<nx; ix++)
		taper[ix]=1.0;
	for(ix=0; ix<range; ix++){
		taper[ix]=(float)(ix+1)/range;
		taper[nx-1-ix]=taper[ix];
	}

	for(i3=0; i3<n3; i3++){
		sf_floatread(dd[0], nt*nx, Fdat);
		for(it=0; it<nt; it++){
			for(ix=0; ix<nx; ix++){
				dd[ix][it] *= taper[ix];
			}
		}
		sf_floatwrite(dd[0], nt*nx, Fnew);
	}

	exit(0);
}
