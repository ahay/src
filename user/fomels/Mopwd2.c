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
    int i, n1, n2, n12;
    float **dat, **ang, **p1, **p2;
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
    p1 = sf_floatalloc2(n1,n2);
    p2 = sf_floatalloc2(n1,n2);

    sf_floatread(dat[0],n12,inp);
    sf_floatread(ang[0],n12,dip);

    for (i=0; i < n12; i++) {
	p1[0][i] = sinf(ang[0][i]);
	p2[0][i] = cosf(ang[0][i]);
    }

    opwd_init(n1,n2);

    opwd_filter(interp,interp,t1,t2,p1,p2,dat,ang);

    for (i=0; i < n12; i++) {
	dat[0][i] -= ang[0][i];
    }

    sf_floatwrite(dat[0],n12,out);

    exit(0);
}
