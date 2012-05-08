/* 2-D omnidirectional plane-wave destruction */
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

#include "opwd.h"

int main(int argc, char *argv[])
{
    char *type;
    sf_tris t1=NULL, t2=NULL;
    interpolate interp;
    int i1, i2, n1, n2, n12;
    float **dat, **ang, **res;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    dip = sf_input("angle"); 
    /* dip angle */

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (NULL == (type = sf_getstring("interp"))) type="lagrange";
    /* interpolation type */

    switch(type[0]) {
	case 'l':
	    interp = lagrange;
	    break;
	case 'b':
	    interp = bspline;
	    t1 = sf_tridiagonal_init(n1);
	    t2 = sf_tridiagonal_init(n2);
	    sf_tridiagonal_const_define(t1,0.75,0.125,true);
	    sf_tridiagonal_const_define(t2,0.75,0.125,true);
	    break;
	default:
	    interp = NULL;
	    sf_error("unknown interp=%s",type);
    }

    dat = sf_floatalloc2(n1,n2);
    ang = sf_floatalloc2(n1,n2);
    res = sf_floatalloc2(n2,n1);

    sf_floatread(dat[0],n12,inp);
    sf_floatread(ang[0],n12,dip);

    filter2(n1,n2,interp,interp,t1,t2,ang,res,dat,ang);

    for (i2=0; i2 < n2; i2++) {
	for (i1=1; i1 < n1-1; i1++) {
	    dat[i2][i1] -= ang[i2][i1];
	}
    }

    sf_floatwrite(dat[0],n12,out);

    exit(0);
}
