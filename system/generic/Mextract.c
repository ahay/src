/* Forward interpolation in 2-D slices. */
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

#include <float.h>
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int id, nk, nd, nm, nt, it, nx, ny, xkey, ykey, interp;
    float *mm, *dd, **xy, *hdr, x0, y0, dx, dy;
    sf_file in, out, head;

    sf_init (argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nx)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&ny)) sf_error("Need n2= in in");
    nt = sf_leftsize(in,2);
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("xkey",&xkey)) xkey=0;
    /* x key number */
    if (!sf_getint("ykey",&ykey)) ykey=1;
    /* y key number */

    /* create coordinates */
    head = sf_input("head");
    if (!sf_histint(head,"n1",&nk)) sf_error("No   n1= in head");
    if (!sf_histint(head,"n2",&nd)) sf_error("Need n2= in head");
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");

    xy  = sf_floatalloc2(2,nd);
    hdr = sf_floatalloc(nk);

    for (id=0; id<nd; id++) {	
	sf_floatread (hdr,nk,head);
	xy[id][0] = hdr[xkey];
	xy[id][1] = hdr[ykey];
    }

    sf_fileclose (head);

    sf_putint(out,"n1",nd);
    sf_putint(out,"n2",1);

    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&y0)) sf_error("No o2= in input");
    
    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* [1,2] interpolation method, 1: nearest neighbor, 2: bi-linear */

    switch (interp) {
	case 1:
	    sf_int2_init (xy, x0,y0,dx,dy,nx,ny, sf_bin_int, 1, nd);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    sf_int2_init (xy, x0,y0,dx,dy,nx,ny, sf_lin_int, 2, nd);
	    sf_warning("Using linear interpolation");
	    break;
	case 3:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (mm,nm,in);

	sf_int2_lop(false,false,nm,nd,mm,dd);

	sf_floatwrite (dd,nd,out);
    }


    exit(0);
}

/* 	$Id: Mextract.c 7107 2011-04-10 02:04:14Z ivlad $	 */
