/* 4-D data binning. 

   if xkey < 0, the first axis indexes traces in a gather like cdpt.
*/
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

#include <limits.h>
#include <string.h>

#include <rsf.h>

#include "segy.h"

int main (int argc, char* argv[])
{
    int id, nk, nd, nt, nx, ny, nz, n2, xkey, ykey, zkey, *hdr, *x, *y, *z;
    int xmin, xmax, ymin, ymax, zmin, zmax, i, j, k, ix, iy, iz, jp=0, kp=0, ***map;
    off_t pos;
    char *buf, *zero, *xk, *yk, *zk, *header;
    sf_file in, out, head, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nd)) sf_error("Need n2= in in");

    header = sf_getstring("head");
    /* header file */
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_INT != sf_gettype(head)) sf_error("Need int header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

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
	/* y key number (if no yk), default is iline */
	ykey = segykey("xline");
    }
    if (NULL != (zk = sf_getstring("zk"))) {
	/* z key name */
	zkey = segykey(zk);
    }  else if (!sf_getint("zkey",&zkey)) {
	/* z key number (if no zk), default is xline */
	zkey = segykey("iline");
    }
    
    if (xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);
    if (zkey < 0 || zkey >= nk) 
	sf_error("zkey=%d is out of the range [0,%d]",zkey,nk-1);

    hdr = sf_intalloc(nk);
    x = sf_intalloc(nd);
    y = sf_intalloc(nd);
    z = sf_intalloc(nd);

    zmin = ymin = xmin = +INT_MAX;
    zmax = ymax = xmax = -INT_MAX;
    for (i=id=0; id<nd; id++) {	
	sf_intread (hdr,nk,head);
	j = hdr[ykey]; 
	if (j < ymin) ymin=j;
	if (j > ymax) ymax=j;
	y[id] = j;
	k = hdr[zkey]; 
	if (k < zmin) zmin=k;
	if (k > zmax) zmax=k;
	z[id] = k;
	if (xkey < 0) { /* index traces in a gather */
	    if (i > 0 && (j != jp || k != kp)) i=0;
	    if (i < xmin) xmin=i;
	    if (i > xmax) xmax=i;
	    x[id] = i;
	    i++;
	    jp = j;
	    kp = k;
	} else {
	    i = hdr[xkey]; 
	    if (i < xmin) xmin=i;
	    if (i > xmax) xmax=i;
	    x[id] = i;
	}
    }

    sf_fileclose (head);

    /* let user overwrite */
    sf_getint ("xmin",&xmin); /* x minimum */
    sf_getint ("xmax",&xmax); /* x maximum */
    sf_getint ("ymin",&ymin); /* y minimum */
    sf_getint ("ymax",&ymax); /* y maximum */
    sf_getint ("zmin",&zmin); /* z minimum */
    sf_getint ("zmax",&zmax); /* z maximum */

    if (xmax < xmin) sf_error ("xmax=%d < xmin=%d",xmax,xmin);
    if (ymax < ymin) sf_error ("ymax=%d < ymin=%d",ymax,ymin);
    if (zmax < zmin) sf_error ("zmax=%d < zmin=%d",zmax,zmin);

    nx = xmax-xmin+1;
    ny = ymax-ymin+1;
    nz = zmax-zmin+1;

    sf_putint (out,"n2",nx);
    sf_putint (out,"n3",ny);
    sf_putint (out,"n4",nz);
    sf_putint (out,"o2",xmin);
    sf_putint (out,"o3",ymin);
    sf_putint (out,"o4",zmin);
    sf_putint (out,"d2",1);
    sf_putint (out,"d3",1);
    sf_putint (out,"d4",1);

    map = sf_intalloc3(nx,ny,nz);

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		map[iz][iy][ix] = -1;
	    }
	}
    }

    for (id = 0; id < nd; id++) {
	ix = x[id];
	iy = y[id];
	iz = z[id];
	if (ix >= xmin && ix <= xmax && 
	    iy >= ymin && iy <= ymax &&
	    iz >= zmin && iz <= zmax) {
	    map[iz-zmin][iy-ymin][ix-xmin] = id;
	}
    }
   
    nt *= sf_esize(in);
	    
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    buf = sf_charalloc(nt);
    zero = sf_charalloc(nt);
    memset(zero,0,nt);

    sf_unpipe(in,(long) nt*nd);

    pos = sf_tell(in);

    /* loop over output */
    for (iz=0; iz < nz; iz++) {
	sf_warning("%d of %d;",iz+1,nz);
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		id = map[iz][iy][ix];
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

    header = sf_getstring("mask");
    /* output mask file */
    if (NULL != header) {
	mask = sf_output(header);
	sf_putint(mask,"n1",nx);
	sf_putint(mask,"n2",ny);
	sf_putint(mask,"n3",nz);
	sf_settype(mask,SF_INT);

	for (iz=0; iz < nz; iz++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    map[iz][iy][ix] = (map[iz][iy][ix] >= 0);
		}
	    }
	}

	sf_intwrite(map[0][0],nx*ny*nz,mask);
    }

    exit(0);
}


