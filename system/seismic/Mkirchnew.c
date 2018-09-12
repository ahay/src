/* Kirchhoff 2-D post-stack time migration and modeling with antialiasing. 

 Antialiasing by reparameterization. */
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

#include <math.h>

#include <rsf.h>

#include "kirchnew.h"

int main(int argc, char* argv[])
{
    int n12, n1, n2, n3, i1, i3, sw;
    bool adj, hd;
    float **data, **modl, *vrms, o1,d1,o2,d2, v0;
    char *test;
    sf_file in, out, vel;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("adj",&adj)) adj=true;
    /* yes: migration, no: modeling */
    if (!sf_getbool("hd",&hd)) hd=true;
    /* if y, apply half-derivative filter */
    if (!sf_getint("sw",&sw)) sw=0;
    /* if > 0, select a branch of the antialiasing operation */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    vrms = sf_floatalloc(n1);

    if (NULL != (test = sf_getstring("velocity"))) { 
	/* velocity file */
	free(test);
	vel = sf_input("velocity");
	sf_floatread(vrms,n1,vel);
	sf_fileclose(vel);
    } else {
	if (!sf_getfloat("v0",&v0)) sf_error("Need velocity= or v0=");
	/* constant velocity (if no velocity=) */
	for (i1=0; i1 < n1; i1++) {
	    vrms[i1] = v0;
	}
    }

    n12 = n1*n2;
    data = sf_floatalloc2(n1,n2);
    modl = sf_floatalloc2(n1,n2);

    kirchnew_init (vrms, o1, d1, d2, n1, n2, sw, true, hd);

    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_floatread (data[0],n12,in);
	} else {
	    sf_floatread (modl[0],n12,in);
	}

	kirchnew_lop (adj,false,n12,n12,modl[0],data[0]);

	if (adj) {
	    sf_floatwrite (modl[0],n12,out);
	} else {
	    sf_floatwrite (data[0],n12,out);
	}
    }


    exit(0);
}
