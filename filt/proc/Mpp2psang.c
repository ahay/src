/* Transform PP angle gathers to PS angle gathers. */
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

#include "int1.h"
#include "prefilter.h"

int main(int argc, char* argv[])
{
    bool inv;
    int nt, na, it, ia, nw, n3, i3;
    float **gather, *trace, *modl, *coord, *gamma; 
    float da, a0, r;
    sf_file in, out, vpvs;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");

    if (!sf_histint(in,"n2",&na)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&da)) sf_error("No d2= in input"); 
    if (!sf_histfloat(in,"o2",&a0)) sf_error("No o2= in input"); 

    prefilter_init (nw,na,2*na);

    n3 = sf_leftsize(in,2);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getint("nw",&nw)) nw=4;
    /* accuracy level */

    gather = sf_floatalloc2(nt,na);
    trace = sf_floatalloc(na);
    coord = sf_floatalloc(na);
    modl =  sf_floatalloc(na);
    gamma = sf_floatalloc(nt);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(gamma,nt,vpvs);
	sf_floatread(gather[0],nt*na,in);

	for (it=0; it < nt; it++) {
	    r = 0.5*(gamma[it]+1./gamma[it]);
	    if (inv) r = 1./r;

	    for (ia=0; ia < na; ia++) {
		coord[ia] = 0.5*acos(r*cos((a0+ia*da)*SF_PI/90.));
		trace[ia] = gather[ia][it];
	    }
	    prefilter_apply (na,trace);
	    int1_init (coord, a0, da, na, sf_spline_int, nw, na);
	    int1_lop (false,false,na,na,trace,modl);

	    for (ia=0; ia < na; ia++) {
		gather[ia][it] = modl[ia];
	    }
	}
	sf_floatwrite(gather[0],nt*na,out);
    }
    
    exit(0);
}

