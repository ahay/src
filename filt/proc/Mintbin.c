/* Data binning. */
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

#include <limits.h>
#include <string.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int id, nk, nd, nt, nx, ny, n2, xkey, ykey, *hdr, *x, *y;
    int xmin, xmax, ymin, ymax, i, ix, iy, **map, esize;
    long pos;
    char *buf, *zero, *xk, *yk;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nd)) sf_error("Need n2= in in");

    head = sf_input("head");

    if (SF_INT != sf_gettype(head)) sf_error("Need int header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    if (NULL != (xk = sf_getstring("xk"))) {
	/* x key name */
	xkey = sf_segykey(xk);
    }  else if (!sf_getint("xkey",&xkey)) {
	/* x key number (if no xk), default is fldr */
	xkey = sf_segykey("fldr");
    }
    if (NULL != (yk = sf_getstring("yk"))) {
	/* y key name */
	ykey = sf_segykey(yk);
    }  else if (!sf_getint("ykey",&ykey)) {
	/* y key number (if no yk), default is tracf */
	ykey = sf_segykey("tracf");
    }
    
    if (xkey < 0 || xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);

    hdr = sf_intalloc(nk);
    x = sf_intalloc(nd);
    y = sf_intalloc(nd);

    ymin = xmin = +INT_MAX;
    ymax = xmax = -INT_MAX;
    for (id=0; id<nd; id++) {	
	sf_intread (hdr,nk,head);
	i = hdr[xkey]; 
	if (i < xmin) xmin=i;
	if (i > xmax) xmax=i;
	x[id] = i;
	i = hdr[ykey]; 
	if (i < ymin) ymin=i;
	if (i > ymax) ymax=i;
	y[id] = i;
    }

    sf_fileclose (head);

    /* let user overwrite */
    sf_getint ("xmin",&xmin);
    sf_getint ("xmax",&xmax);
    sf_getint ("ymin",&ymin);
    /* Grid dimensions */
    sf_getint ("ymax",&ymax);

    if (xmax < xmin) sf_error ("xmax=%d < xmin=%d",xmax,xmin);
    if (ymax < ymin) sf_error ("ymax=%d < ymin=%d",xmax,xmin);

    nx = xmax-xmin+1;
    ny = ymax-ymin+1;

    sf_putint (out,"n2",nx);
    sf_putint (out,"n3",ny);
    sf_putint (out,"o2",xmin);
    sf_putint (out,"o3",ymin);
    sf_putint (out,"d2",1);
    sf_putint (out,"d3",1);

    map = sf_intalloc2(nx,ny);

    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    map[iy][ix] = -1;
	}
    }

    for (id = 0; id < nd; id++) {
	ix = x[id];
	iy = y[id];
	if (ix >= xmin && ix <= xmax && 
	    iy >= ymin && iy <= ymax) {
	    map[iy-ymin][ix-xmin] = id;
	}
    }

    if (!sf_histint(in,"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("wrong esize=%d",esize);
    }
    nt *= esize;
	    
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    buf = sf_charalloc(nt);
    zero = sf_charalloc(nt);
    memset(zero,0,nt);

    sf_unpipe(in,(long) nt*nd);

    pos = sf_tell(in);

    /* loop over output */
    for (iy=0; iy < ny; iy++) {
	sf_warning("%d of %d",iy+1,ny);
	for (ix=0; ix < nx; ix++) {
	    id = map[iy][ix];
	    if (id < 0) {
		sf_charwrite (zero,nt,out);
	    } else {
		sf_seek(in,pos + (long) id*nt,SEEK_SET);		
		sf_charread (buf,nt,in);
		sf_charwrite (buf,nt,out);
	    }
	}
    }

    sf_close();
    exit(0);
}

/* 	$Id$	 */
