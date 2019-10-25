/* Data regularization in 3-D using plane-wave destruction. */
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
#include "int2.h"

int main(int argc, char* argv[])
{
    int niter, nw, n1, n2, n123, nj1, nj2, nk, xkey, ykey, nx, ny;
    int n3, i3, interp, nt, id, nd;
    float *mm, *dd, *pp, *qq, **xy, x0, dx, y0, dy, eps, *hdr;
    char *header;
    bool verb, drift;
    sf_file in, out, dip, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    /* irregular data */
    if (!sf_histint(in,"n1",&nt)) nt=1;
    if (!sf_histint(in,"n2",&nd)) nd=1;
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    n3 = sf_leftsize(in,2);

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);
    
    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    if (!sf_getint("xkey",&xkey)) sf_error("Need xkey=");
    /* x key number */
    if (!sf_getint("ykey",&ykey)) sf_error("Need ykey=");
    /* y key number */

    hdr = sf_floatalloc(nk);
    for (id=0; id<nd; id++) {	
	sf_floatread (hdr,nk,head);
	xy[id][0] = hdr[xkey];
	xy[id][1] = hdr[ykey]; 
    }
    sf_fileclose (head);


    /* regular data */
    if (!sf_histint(dip,"n1",&n1) || n1 != nt) sf_error("Need n1=%d in dip",nt);
    if (!sf_histint(dip,"n2",&nx)) sf_error("Need n2= in dip");
    if (!sf_histint(dip,"n3",&ny)) sf_error("Need n3= in dip");
    if (SF_FLOAT != sf_gettype(dip)) sf_error("Need float data in dip");

    if (!sf_histfloat(dip,"o2",&x0)) sf_error("Need o2= in dip");
    if (!sf_histfloat(dip,"d2",&dx)) sf_error("Need d2= in dip");
    if (!sf_histfloat(dip,"o3",&y0)) sf_error("Need o3= in dip");
    if (!sf_histfloat(dip,"d3",&dy)) sf_error("Need d3= in dip");

    sf_putint(out,"n2",nx);
    sf_putfloat(out,"o2",x0);
    sf_putfloat(out,"d2",dx);
    sf_putint(out,"n3",ny);
    sf_putfloat(out,"o3",y0);
    sf_putfloat(out,"d3",dy);
    sf_putint(out,"n4",n3);

    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    int2_init (xy, x0,y0,dx,dy,nx,ny, nt, sf_spline_int, interp, nd);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("nj1",&nj1)) nj1=1;
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    n123 = nt*nx*ny;
    
    pp = sf_floatalloc(n123);
    qq = sf_floatalloc(n123);

    mm = sf_floatalloc(n123);
    dd = sf_floatalloc(nt*nd);

    allpass3_init(allpass_init(nw, nj1, nt,nx,ny, drift, pp),
		  allpass_init(nw, nj2, nt,nx,ny, drift, qq));

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(pp,n123,dip);
	sf_floatread(qq,n123,dip);

	sf_floatread(dd,nt*nd,in);

	sf_solver_reg(int2_lop, sf_cgstep, allpass3_lop, 2*n123, n123, nt*nd, mm, dd, niter, eps, "verb", verb, "end");
	sf_cgstep_close();

	sf_floatwrite (mm,n123,out);
    }


    exit(0);
}
