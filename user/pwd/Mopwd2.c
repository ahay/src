/* 2-D omnidirectional plane wave destruction. */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "opwd2.h"

int main (int argc, char *argv[])
{
    int ir, nr, n1,n2, m1, m2, n12, nw, i3;
    float *u1, *u2, *p1, *p2;
    sf_file in, out, ang;
    omni2 ap;

    sf_init(argc,argv);
    in = sf_input ("in");
    ang = sf_input ("angle");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(ang)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n12 = n1*n2;
    nr = sf_leftsize(in,2);

    if (!sf_histint(ang,"n1",&m1) || m1 != n1) 
	sf_error("Need n1=%d in angle",n1);
    if (1 != n2 && (!sf_histint(ang,"n2",&m2) || m2 != n2)) 
	sf_error("Need n2=%d in angle",n2);
    
    if (!sf_getint("order",&nw)) nw=1;
    /* accuracy */

    u1 = sf_floatalloc(n12);
    u2 = sf_floatalloc(n12);
    p1 = sf_floatalloc(n12);
    p2 = sf_floatalloc(n12);

    for (ir=0; ir < nr; ir++) {
	/* read dip */
	sf_floatread(u1,n12,ang);
	for (i3=0; i3 < n12; i3++) {
	    p1[i3] = sinf(u1[i3]);
	    p2[i3] = cosf(u1[i3]);
	}
	
	/* read data */
	sf_floatread(u1,n12,in);
		
	ap = opwd2_init (nw,n1,n2,p1,p2);
		
	/* apply */
	opwd21(false, ap, u1, u2);
		
	/* write out */
	sf_floatwrite(u2,n12,out);
    }
	        
    exit (0);
}
