/* Smoothing by spraying */
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

#include "smspray.h"

int main(int argc, char* argv[])
{
    bool adj;
    char *type;
    int i2, n1, n2, ns;
    float *trace, *smooth;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    if (!sf_getbool("adj",&adj)) adj=false; 
    /* adjoint flag */

    if (!sf_getint("ns",&ns)) ns=0;
    /* smoothing radius */

    if (NULL==(type=sf_getstring("type"))) type="triangle";
    /* weight type (triangle, gauss) */

    smspray_init(n1,ns,type[0]);

    trace = sf_floatalloc(n1);
    smooth = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,inp);

	if (adj) {
	    smspray_lop(true,false,n1,n1,smooth,trace);
	} else {
	    smspray_lop(false,false,n1,n1,trace,smooth);
	}
    
	sf_floatwrite(smooth,n1,out);
    }

    exit(0);
}
