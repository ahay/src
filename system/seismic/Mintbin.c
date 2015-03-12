/* Data binning by trace sorting. 

If inv=n, the input is 2-D (n1 x ntr). The output is 3-D (n1 x n2 x n3), n2 and
n3 correspond to two selected keys from the header file. 

If inv=y, the input is 3-D, and the output is 2-D.

if xkey < 0, the first axis indexes traces in a gather like cdpt.
*/
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

#include "segy.h"

int main (int argc, char* argv[])
{
    bool inv;
    int id, nk, nd, nt, nx, ny, n2, xkey, ykey, *hdr, *x, *y;
    int xmin, xmax, ymin, ymax, i, ix, iy, **map, nxy, j, jp=0;
    off_t pos;
    char *buf, *zero, *xk, *yk, *header;
    sf_file in, out, head, mask, mapf;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false; 
    /* inversion flag */

    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1= in in");

    header = sf_getstring("head");
    /* header file */
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_INT != sf_gettype(head)) sf_error("Need int header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    n2 = sf_leftsize(head,1);

    if (inv) {
	nd = n2;
    } else {
	nd = sf_leftsize(in,1);
	if (n2 != nd) sf_error("Wrong number of traces in head");
    }

    segy_init(nk,head);

    if (NULL != (xk = sf_getstring("xk"))) {
	/* x key name */
	xkey = segykey(xk);
    }  else if (!sf_getint("xkey",&xkey)) {
	/* x key number (if no xk), default is fldr */
	xkey = segykey("fldr");
    }
    if (NULL != (yk = sf_getstring("yk"))) {
	/* y key name */
	ykey = segykey(yk);
    }  else if (!sf_getint("ykey",&ykey)) {
	/* y key number (if no yk), default is tracf */
	ykey = segykey("tracf");
    }
    
    if (xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);

    hdr = sf_intalloc(nk);
    x = sf_intalloc(nd);
    y = sf_intalloc(nd);

    for (i=id=0; id<nd; id++) {	
	sf_intread (hdr,nk,head);
	j = hdr[ykey];
	y[id] = j; 
	if (xkey < 0) { /* index traces in a gather */
	    if (i > 0 && j != jp) i=0;
	    x[id] = i;
	    i++;
	    jp = j;
	} else {
	    i = hdr[xkey];
	    x[id] = i;
	}
    }

    sf_fileclose (head);

    if (inv) {	
	if (!sf_histint (in,"n2",&nx))   sf_error("No n2= in input");
	if (!sf_histint (in,"n3",&ny))   sf_error("No n3= in input");
	if (!sf_histint (in,"o2",&xmin)) sf_error("No o2= in input");
	if (!sf_histint (in,"o3",&ymin)) sf_error("No o3= in input");

	xmax = xmin+nx-1;
	ymax = ymin+ny-1;

	sf_putint (out,"n2",nd);
	sf_putint (out,"n3",1);
    } else {
	ymin = xmin = +INT_MAX;
	ymax = xmax = -INT_MAX;
	for (id=0; id<nd; id++) {	
	    i = x[id]; 
	    if (i < xmin) xmin=i;
	    if (i > xmax) xmax=i;
	    i = y[id];
	    if (i < ymin) ymin=i;
	    if (i > ymax) ymax=i;
	}
	
	/* let user overwrite */
	sf_getint ("xmin",&xmin); /* x minimum */
	sf_getint ("xmax",&xmax); /* x maximum */
	sf_getint ("ymin",&ymin); /* y minimum */
	sf_getint ("ymax",&ymax); /* y maximum */
	
	if (xmax < xmin) sf_error ("xmax=%d < xmin=%d",xmax,xmin);
	if (ymax < ymin) sf_error ("ymax=%d < ymin=%d",ymax,ymin);

	nx = xmax-xmin+1;
	ny = ymax-ymin+1;
	
	sf_putint (out,"n2",nx);
	sf_putint (out,"n3",ny);
	sf_putint (out,"o2",xmin);
	sf_putint (out,"o3",ymin);
	sf_putint (out,"d2",1);
	sf_putint (out,"d3",1);
    }

    nxy = nx*ny;

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
   
    nt *= sf_esize(in);
	    
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    buf = sf_charalloc(nt);

    zero = sf_charalloc(nt);
    memset(zero,0,nt);

    if (inv) {
	sf_unpipe(in,(off_t) nt*nxy); 
    } else { 
	sf_unpipe(in,(off_t) nt*nd);    
    }

    pos = sf_tell(in);

    /* loop over output */
    if (inv) {
	for (id = 0; id < nd; id++) {
	    ix = x[id]-xmin;
	    iy = y[id]-ymin;
	    if (ix < 0 || ix >= nx || 
		iy < 0 || iy >= ny) {
		sf_charwrite (zero,nt,out);
	    } else {
		sf_seek(in,pos + (off_t) (ix+iy*nx)*nt,SEEK_SET);		
		sf_charread (buf,nt,in);
		sf_charwrite (buf,nt,out);
	    }
	}
    } else {
	for (iy=0; iy < ny; iy++) {
	    sf_warning("%d of %d;",iy+1,ny);
	    for (ix=0; ix < nx; ix++) {
		id = map[iy][ix];
		
		if (id < 0) {
		    sf_charwrite (zero,nt,out);
		} else {
		    sf_seek(in,pos + (off_t) id*nt,SEEK_SET);		
		    sf_charread (buf,nt,in);
		    sf_charwrite (buf,nt,out);
		}
	    }
	}
    }
    sf_warning(".");

    header = sf_getstring("map");
    /* output map file */
    if (NULL != header) {
	mapf = sf_output(header);
	sf_putint(mapf,"n1",nx);
	sf_putint(mapf,"n2",ny);
	sf_putint(mapf,"n3",1);
	sf_settype(mapf,SF_INT);

	sf_intwrite(map[0],nxy,mapf);
    }

    header = sf_getstring("mask");
    /* output mask file */
    if (NULL != header) {
	mask = sf_output(header);
	sf_putint(mask,"n1",nx);
	sf_putint(mask,"n2",ny);
	sf_putint (mask,"o1",xmin);
	sf_putint (mask,"o2",ymin);
	sf_putint (mask,"d1",1);
	sf_putint (mask,"d2",1);
	sf_putstring (mask,"unit1","");
	sf_putstring (mask,"unit2","");
	sf_putstring (mask,"label1",(NULL==xk)? "":xk);
	sf_putstring (mask,"label2",(NULL==yk)? "":yk);
	sf_settype(mask,SF_INT);

	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		map[iy][ix] = (map[iy][ix] >= 0);
	    }
	}

	sf_intwrite(map[0],nxy,mask);
    }


    exit(0);
}

/* 	$Id$	 */
