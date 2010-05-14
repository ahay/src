/* Test linear matching operator. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "matching.h"       

int main (int argc, char **argv)
{
    bool verb;                /* verbosity flag */
    bool adj;
    int nt;                   /* number of time samples */
    int nx;                   /* number of offsets */ 
    int nf1, nf2, nf;         /* number of filter size */
    int n3, i3;
    float *x, *y, *filt;      /* input and output */
    sf_file in, out, fil;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, perform adjoint operation */

    if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity flag */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histint(fil,"n1",&nf1)) sf_error("No nf1= in filt");
    if (!sf_histint(fil,"n2",&nf2)) sf_error("No nf2= in filt");
    if (!sf_histint(fil,"n3",&nf)) sf_error("No nf3= in filt");

    if (nt != nf1 || nx != nf2 ) sf_error("Nee n1==nf1 && n2==nf2");

    x = sf_floatalloc (nt*nx);
    y = sf_floatalloc (nt*nx);
    filt = sf_floatalloc (nt*nx*nf);

    for (i3=0; i3 < n3; i3++) {
	if(verb) sf_warning("i=%d of %d",i3+1,n3);

	sf_floatread(filt,nt*nx*nf,fil);

	matching_init(filt,nt,nx,nf);

	if (adj) {
	    sf_floatread(y,nt*nx,in);
	} else { /* modeling */
	    sf_floatread(x,nt*nx,in);
	}
    
	pmatching_lop (adj, false, nt*nx, nt*nx, x, y);
	
	if (adj) {
	    sf_floatwrite(x,nt*nx,out);
	} else { /* modeling */
	    sf_floatwrite(y,nt*nx,out);
	}
    }
    
    exit (0);
}

/* 	$Id$	 */
