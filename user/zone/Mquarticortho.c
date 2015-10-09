/* Interval quartic coefficients estimation
*/
/*
 Copyright (C) 2009 University of Texas at Austin
 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Main program------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	int nz,nx,ny,nlayer;
	float dx,ox,dy,oy,dz,oz;
	float **eff,*twtime, **inverted, ***mod;

	sf_file model, inv, ef, efft;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input-----------------------------------------------------------------------------------------*/
	model = sf_input("in"); /* reflector model*/
	if (!sf_histint(model,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histint(model,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histint(model,"n2",&ny)) sf_error("No n2= in input");
	
	if (!sf_histfloat(model,"o1",&oz)) oz=0.;
	if (!sf_histfloat(model,"o2",&ox)) ox=0.;
	if (!sf_histfloat(model,"o3",&oy)) oy=0.;
	
	if (!sf_histfloat(model,"d1",&dz)) dz=1.;
	if (!sf_histfloat(model,"d2",&dx)) dx=1.;
	if (!sf_histfloat(model,"d3",&dy)) dy=1.;
	
	ef = sf_input("eff"); /* effective parameters at the bottom of each layer*/
	efft = sf_input("twtime"); /* effective two-way traveltime at the bottom of each layer*/
	
	if (!sf_histint(thickness,"n1",&nlayer)) sf_error("No n1= in thickness");
	
	/* Allocate space-------------------------------------------------------------------------------------*/
	mod = sf_floatalloc3(nz,nx,ny); /* model*/
	eff = sf_floatalloc2(8,nlayer);
	inverted = sf_floatalloc2(8,nlayer);
	twtime = sf_floatalloc(nlayer);
	
	/* Read input-----------------------------------------------------------------------------------------*/
	sf_floatread (eff[0],8*nlayer,ef);
	sf_floatread (twtime,nlayer,efft);
	
	/* Set Output------------------------------------------------------------------------------------------*/
	inv = sf_output("out");
	sf_putint(inv,"n1",8);
	sf_putint(inv,"d1",1);
	sf_putint(inv,"o1",0);
	
	sf_putint(inv,"n2",nlayer);
	sf_putint(inv,"d2",1);
	sf_putint(inv,"o2",0);
	
	/* Inversion-------------------------------------------------------------------------------------------*/
	/* Step 1 : Generalized Dix----------------------------------------------------------------------------*/
	
	gendix(inverted, eff, twtime, nlayer);
	
	/* Step 2 : Quartic Inversion---------------------------------------------------------------------------*/
	
	quartic(inverted, eff, twtime, nlayer);
	
	/* Write output*/
	sf_floatwrite(inverted[0],8*nlayer,inv);
	
	exit(0);
}
