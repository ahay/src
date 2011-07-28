/* Simple operations with real sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "rsin.h"

int main(int argc, char* argv[])
{
    int n1, i2, n2;
    bool adj;
    const char *op;
    float *c0, *x, *y;
    sf_operator oper=NULL;
    sf_file in, out, root;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    root = sf_input("root");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (sf_filesize(root) != n2) sf_error("Wrong dimensions in root");
    c0 = sf_floatalloc(n2);
    sf_floatread(c0,n2,root);
    sf_fileclose(root);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (NULL == (op = sf_getstring("oper"))) op="destruct";
    /* operation to perform */

    switch (op[0]) {
	case 'd':
	    oper = rsin_destruct;
	    break;
	default:
	    sf_error("Unknown operator \"%s\"",op);
    }

    x = sf_floatalloc(n1);
    y = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(x,n1,in);

	rsin_init(c0[i2]);

	if (adj) {
	    oper(true,false,n1,n1,y,x);
	} else {
	    oper(false,false,n1,n1,x,y);
	}

	sf_floatwrite(y,n1,out);
    }


    exit(0);
}
    
