/* ENO interpolation in 2-D slices. */
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

int main (int argc, char* argv[])
{
    int id, nk, nd, nm, nt, it, nx, ny, xkey, ykey, interp, i, j;
    float *mm, *dd, **xy, *hdr, x0, y0, dx, dy, x, y, dt, t0, f1[2];
    char *xk, *yk;
    sf_eno2 ent;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    head = sf_input("head");

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&nd)) sf_error("Need n1= in in");

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&ny)) sf_error("Need n2= in in");
    if (!sf_histint(in,"n3",&nt)) sf_error("Need n3= in in");

    if (!sf_histfloat(in,"d1",&dx)) dx=1.;
    if (!sf_histfloat(in,"d2",&dy)) dy=1.;
    if (!sf_histfloat(in,"o1",&x0)) x0=0.;
    if (!sf_histfloat(in,"o2",&y0)) y0=0.;

    sf_putint(out,"n1",nd);    
    sf_putint(out,"n2",nt);
    sf_putint(out,"n3",1);

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

    hdr = sf_floatalloc(nk);
    for (id=0; id<nd; id++) {	
	sf_floatread (hdr,nk,head);
	xy[id][0] = hdr[xkey];
	xy[id][1] = hdr[ykey];
    }
    sf_fileclose (head);

    if (sf_histfloat(in,"o3",&t0)) sf_putfloat(out,"o2",t0);
    if (sf_histfloat(in,"d3",&dt)) sf_putfloat(out,"d2",dt);

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation order */

    ent = sf_eno2_init (interp, nx, ny);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (mm,nm,in);
	sf_eno2_set1 (ent,mm);
	for (id=0; id<nd; id++) {
	    x = (xy[id][0]-x0)/dx; i=x; x-= i;
	    y = (xy[id][1]-y0)/dy; j=y; y-= j;
	    sf_eno2_apply(ent, i, j, x, y, dd+id, f1, FUNC);
	}
	sf_floatwrite (dd,nd,out);
    }

    exit(0);
}

/* 	$Id: Mbin.c 847 2004-10-27 20:33:16Z fomels $	 */
