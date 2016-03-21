/* Merging legacy and hires datasets */
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

#include "nsmooth1.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n1, n2, n12;
    float *legacy, *hires, *merge, *wght, **nr;
    sf_file in, out, weight, rect;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    rect = sf_input("rect");
    weight = sf_input("weight");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(rect)) sf_error("Need float rect");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    legacy = sf_floatalloc(n12);
    hires = sf_floatalloc(n12);
    merge = sf_floatalloc(n12);
    wght = sf_floatalloc(n12);
    nr = sf_floatalloc2(n1,n2);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    sf_floatread(nr[0],n12,rect);
    sf_floatread(wght,n12,weight);
    
    nsmooth1_init(n1,n2,nr);
    sf_weight_init(wght);

    if (adj) {
	sf_floatread(hires,n12,in);
	sf_floatread(legacy,n12,in);

	sf_weight_lop(true,false,n12,n12,merge,hires);
	nsmooth1_lop(true,true,n12,n12,merge,legacy);

	sf_floatwrite(merge,n12,out);
    } else {
	sf_floatread(merge,n12,in);

	sf_weight_lop(false,false,n12,n12,merge,hires);
	nsmooth1_lop(false,false,n12,n12,merge,legacy);

	sf_floatwrite(hires,n12,out);
	sf_floatwrite(legacy,n12,out);
    }

    exit(0);
}
	
	
