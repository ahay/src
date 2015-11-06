/* One side of plane wave destruction. */
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

#include "allp3.h"

int main (int argc, char *argv[])
{
    bool left;
    int ir, nr, n1,n2,n3, m1, m2, m3, n12, n123, nw, nj1, nj2, i3;
    float *u1, *u2, *p;
    sf_file in, out, dip;
    off_t pos=0;
    allpass ap;

    sf_init(argc,argv);
    in = sf_input ("in");
    dip = sf_input ("dip");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(dip)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n12 = n1*n2;
    n123 = n12*n3;
    nr = sf_leftsize(in,3);

    if (!sf_getbool("left",&left)) left=true;
    /* if using left or right side of PWD */

    if (!sf_histint(dip,"n1",&m1) || m1 != n1) 
	sf_error("Need n1=%d in dip",n1);
    if (1 != n2 && (!sf_histint(dip,"n2",&m2) || m2 != n2)) 
	sf_error("Need n2=%d in dip",n2);
    if (1 != n3 && (!sf_histint(dip,"n3",&m3) || m3 != n3)) 
	sf_error("Need n3=%d in dip",n3);

    if (!sf_getint("order",&nw)) nw=1;
    /* accuracy */
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* in-line aliasing */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* cross-line aliasing */

    for (ir=0; ir < nr; ir++) {
	u1 = sf_floatalloc(n12);
	u2 = sf_floatalloc(n12);
	p  = sf_floatalloc(n12);
	    
	for (i3=0; i3 < n3; i3++) {
	    /* read data */
	    sf_floatread(u1,n12,in);
	    
	    /* read t-x dip */
	    sf_floatread(p,n12,dip);
		
	    ap = allpass_init (nw,nj1,n1,n2,1,p);
		
	    /* apply */
	    if (left) {
		left1(false, false, ap, u1, u2);
	    } else {
		right1(false, false, ap, u1, u2);
	    }
		
	    /* write t-x destruction */
	    sf_floatwrite(u2,n12,out);
	}
	
	free(u1);
	free(u2);
	free(p);
	    
    }
    
    exit (0);
}

