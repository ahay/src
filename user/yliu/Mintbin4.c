/* 5-D data binning. */
/*
  Copyright (C) 2013 Jilin University
  
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
    int id, nk, nd, nt, nx, ny, ns, nr, n2, dx, dy, ds, dr, tmpd;
    int xkey, ykey, skey, rkey, *hdr, *x, *y, *s, *r;
    int xmin, xmax, ymin, ymax, smin, smax, rmin, rmax;
    int i, ix, iy, is, ir, ****map;
    off_t pos;
    char *buf, *zero, *xk, *yk, *sk, *rk, *header;
    sf_file in, out, head, fold;

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
	/* y key number (if no yk), default is tracf */
	ykey = segykey("tracf");
    }
    if (NULL != (sk = sf_getstring("sk"))) {
	/* s key name */
	skey = segykey(sk);
    }  else if (!sf_getint("skey",&skey)) {
	/* s key number (if no sk), default is swdep */
	skey = segykey("swdep");
    }
    if (NULL != (rk = sf_getstring("rk"))) {
	/* r key name */
	rkey = segykey(rk);
    }  else if (!sf_getint("rkey",&rkey)) {
	/* r key number (if no rk), default is gwdep */
	rkey = segykey("gwdep");
    }
    
    if (xkey < 0 || xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);
    if (skey < 0 || skey >= nk) 
	sf_error("skey=%d is out of the range [0,%d]",skey,nk-1);
    if (rkey < 0 || rkey >= nk) 
	sf_error("rkey=%d is out of the range [0,%d]",rkey,nk-1);

    hdr = sf_intalloc(nk);
    x = sf_intalloc(nd);
    y = sf_intalloc(nd);
    s = sf_intalloc(nd);
    r = sf_intalloc(nd);


    rmin = smin = ymin = xmin = +INT_MAX;
    rmax = smax = ymax = xmax = -INT_MAX;
    dr = ds = dy = dx = tmpd = +INT_MAX;
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
	i = hdr[skey]; 
	if (i < smin) smin=i;
	if (i > smax) smax=i;
	s[id] = i;
	i = hdr[rkey]; 
	if (i < rmin) rmin=i;
	if (i > rmax) rmax=i;
	r[id] = i;
    }

    sf_fileclose (head);

    /* let user overwrite */
    sf_getint ("xmin",&xmin); /* x minimum */
    sf_getint ("xmax",&xmax); /* x maximum */
    sf_getint ("ymin",&ymin); /* y minimum */
    sf_getint ("ymax",&ymax); /* y maximum */
    sf_getint ("smin",&smin); /* s minimum */
    sf_getint ("smax",&smax); /* s maximum */
    sf_getint ("rmin",&rmin); /* r minimum */
    sf_getint ("rmax",&rmax); /* r maximum */

    if (xmax < xmin) sf_error ("xmax=%d < xmin=%d",xmax,xmin);
    if (ymax < ymin) sf_error ("ymax=%d < ymin=%d",ymax,ymin);
    if (smax < smin) sf_error ("smax=%d < smin=%d",smax,smin);
    if (rmax < rmin) sf_error ("rmax=%d < rmin=%d",rmax,rmin);

    if (!sf_getint("dx",&dx)) dx=1;
    if (!sf_getint("dy",&dy)) dy=1;
    if (!sf_getint("ds",&ds)) ds=1;
    if (!sf_getint("dr",&dr)) dr=1;

    nx = (xmax-xmin)/dx+1;
    ny = (ymax-ymin)/dy+1;
    ns = (smax-smin)/ds+1;
    nr = (rmax-rmin)/dr+1;

    sf_putint (out,"n2",nx);
    sf_putint (out,"n3",ny);
    sf_putint (out,"n4",ns);
    sf_putint (out,"n5",nr);

    sf_putint (out,"o2",xmin);
    sf_putint (out,"o3",ymin);
    sf_putint (out,"o4",smin);
    sf_putint (out,"o5",rmin);

    sf_putint (out,"d2",dx);
    sf_putint (out,"d3",dy);
    sf_putint (out,"d4",ds);
    sf_putint (out,"d5",dr);

    map = sf_intalloc4(nx,ny,ns,nr);

    for (ir=0; ir < nr; ir++) {
	for (is=0; is < ns; is++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    map[ir][is][iy][ix] = -1;
		}
	    }
	}
    }

    for (id = 0; id < nd; id++) {
	ix = x[id];
	iy = y[id];
	is = s[id];
	ir = r[id];
	if (ix >= xmin && ix <= xmax && 
	    iy >= ymin && iy <= ymax &&
	    is >= smin && is <= smax &&
	    ir >= rmin && ir <= rmax) {
	    map[(ir-rmin)/dr][(is-smin)/ds][(iy-ymin)/dy][(ix-xmin)/dx] = id;
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
    for (ir=0; ir < nr; ir++) {
	for (is=0; is < ns; is++) {
	    sf_warning("%d of %d, %d of %d;",ir+1,nr,is+1,ns);
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    id = map[ir][is][iy][ix];
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
    }
    sf_warning(".");

    header = sf_getstring("fold");
    /* output fold file */
    if (NULL != header) {
	fold = sf_output(header);
	sf_putint(fold,"n1",nx);
	sf_putint(fold,"n2",ny);
	sf_putint(fold,"n3",ns);
	sf_putint(fold,"n4",nr);
	sf_settype(fold,SF_INT);

	for (ir=0; ir < nr; ir++) {
	    for (is=0; is < ns; is++) {
		for (iy=0; iy < ny; iy++) {
		    for (ix=0; ix < nx; ix++) {
			map[ir][is][iy][ix] = 0;
		    }
		}
	    }
	}
	
	for (id = 0; id < nd; id++) {
	    ix = x[id];
	    iy = y[id];
	    is = s[id];
	    ir = r[id];
	    if (ix >= xmin && ix <= xmax && 
		iy >= ymin && iy <= ymax &&
		is >= smin && is <= smax &&
		ir >= rmin && ir <= rmax) {
		map[(ir-rmin)/dr][(is-smin)/ds][(iy-ymin)/dy][(ix-xmin)/dx] ++;
	    }
	}

	sf_intwrite(map[0][0][0],nx*ny*ns*nr,fold);
    }

    exit(0);
}


