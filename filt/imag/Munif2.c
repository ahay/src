/* Generate layered velocity model from specified interfaces. 

In each layer, velocity is a linear function of position.

Inspired by SU's unif2.
*/
/*
Copyright (C) 2004 University of Texas at Austin

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

int main(int argc, char **argv)
{
    int n1, n2, ninf, i1, i2, i;
    float o1, d1, o2, d2, *v0, *dvdx, *dvdz, *x0, *z0, *trace, **inter, x, z;
    sf_file model, surface;
    
    sf_init(argc, argv);
    surface = sf_input("in");
    model = sf_output("out");

    if (SF_FLOAT != sf_gettype(surface)) sf_error("Need float input");

    if (!sf_histint(surface,"n1",&n2))   sf_error("No n1= in input");
    if (!sf_histfloat(surface,"d1",&d2)) sf_error("No d1= in input");
    if (!sf_histfloat(surface,"o1",&o2)) o2=0.;

    sf_putint(model,"n2",n2);
    sf_putfloat(model,"d2",d2);
    sf_putfloat(model,"o2",o2);

    if (!sf_histint(surface,"n2",&ninf)) ninf=1; 

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* Number of samples on the depth axis */
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* Sampling of the depth axis */
    if (!sf_getfloat("o1",&o1)) o1=0.;
    /* Origin of the depth axis */
	 	
    sf_putint(model,"n1",n1);
    sf_putfloat(model,"d1",d1);
    sf_putfloat(model,"o1",o1);

    inter = sf_floatalloc2(n2,ninf);
    sf_floatread(inter[0],n2*ninf,surface);

    ninf++; /* more layers than interfaces */
    v0 = sf_floatalloc(ninf);
    x0 = sf_floatalloc(ninf);
    z0 = sf_floatalloc(ninf);
    dvdx = sf_floatalloc(ninf);
    dvdz = sf_floatalloc(ninf);

    /* Input layer velocities and velocity derivatives */
    if (!sf_getfloats("x0",x0,ninf)) 
	for(i=0;i< ninf;i++) x0[i] = 0.;
    if (!sf_getfloats("z0",z0,ninf))
	for(i=0;i< ninf;i++) z0[i] = 0.;
    if (!sf_getfloats("v00",v0,ninf))
	for(i=0;i< ninf;i++) v0[i] = 1500.+ 500*i;
    if (!sf_getfloats("dvdx",dvdx,ninf)) 
	for(i=0;i< ninf;i++) dvdx[i] = 0.;
    if (!sf_getfloats("dvdz",dvdz,ninf)) 
	for(i=0;i< ninf;i++) dvdz[i] = 0.;

    trace = sf_floatalloc(n1);

    /* compute linear velocity */
    for(i2=0; i2 < n2; i2++) { 
	x = o2+i2*d2;
	for(i1=0; i1 < n1; i1++) {
	    z = o1 + i1*d1;
	    for (i=0; i < ninf-1; i++) {
		if (z < inter[i][i2]) {
		    trace[i1] = v0[i] + (x-x0[i])*dvdx[i] + (z-z0[i])*dvdz[i];
		    break;
		}
	    }
	    if (i == ninf-1) /* bottom layer */
		trace[i1] = v0[i] + (x-x0[i])*dvdx[i] + (z-z0[i])*dvdz[i];
	}
	sf_floatwrite(trace,n1,model);
    }

    exit(0);
}

/* 	$Id: Munif2.c,v 1.8 2004/07/02 18:08:39 fomels Exp $	 */
