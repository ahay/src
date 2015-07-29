/* Data binning in 2-D slices. 

December 2014 program of the month:
http://ahay.org/rsflog/index.php?/archives/419-Program-of-the-month-sfbin.html
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

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "medbin.h"

int main (int argc, char* argv[])
{
    bool norm;
    int id, nk, nd, im, nm, nt, it, nx, ny, n2, xkey, ykey, interp, i4, n4;
    float *mm, *count, *dd, ***xy, *hdr;
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, clip;
    char *header;
    sf_file in, out, head, fold;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nt)) sf_error("Need n2= in in");
    n4 = sf_leftsize(in,2);

    /* create coordinates */
    xy = sf_floatalloc3(2,nd,n4);

    header = sf_getstring("head");
    /* header file */
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");    
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    if (!sf_getint("xkey",&xkey)) xkey=0;
    /* x key number */

    if (!sf_getint("ykey",&ykey)) ykey=1;
    /* y key number */

    if (xkey < 0 || xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);

    hdr = sf_floatalloc(nk);

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    for (i4=0; i4 < n4; i4++) {
	for (id=0; id<nd; id++) {	
	    sf_floatread (hdr,nk,head);
	    f = hdr[xkey]; 
	    if (f < xmin) xmin=f;
	    if (f > xmax) xmax=f;
	    xy[i4][id][0] = f;
	    f = hdr[ykey]; 
	    if (f < ymin) ymin=f;
	    if (f > ymax) ymax=f;
	    xy[i4][id][1] = f;
	}
    }
    sf_fileclose (head);

    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(out,"o3",t0);
    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(out,"d3",dt);

    /* let user overwrite */
    sf_getfloat ("xmax",&xmax); /* x maximum */
    sf_getfloat ("xmin",&xmin); /* x minimum */
    sf_getfloat ("ymax",&ymax); /* y maximum */
    sf_getfloat ("ymin",&ymin); /* y minimum */
     
    if (xmax < xmin) sf_error ("xmax=%f < xmin=%f",xmax,xmin);
    if (ymax < ymin) sf_error ("ymax=%f < ymin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; /* x origin */
    if (!sf_getfloat("y0",&y0)) y0=ymin; /* y origin */

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    /* create model */
    if (!sf_getint ("nx",&nx)) nx = (int) (xmax - xmin + 1.);
    /* Number of bins in x */
    if (!sf_getint ("ny",&ny)) ny = (int) (ymax - ymin + 1.);
    /* Number of bins in y */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_shiftdim(in, out, 2);

    if (!sf_getfloat("dx",&dx)) {
	/* bin size in x */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }

    if (!sf_getfloat("dy",&dy)) {
	/* bin size in y */
	if (1 >= nx) {
	    dy = dx;
	} else {
	    dy = (ymax-ymin)/(ny-1);
	}
    }

    sf_putfloat (out,"d1",dx);
    sf_putfloat (out,"d2",dy);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=1;
    /* [0,1,2] interpolation method;
       0: median, 1: nearest neighbor, 2: bi-linear,  */

    if (!sf_getbool("norm",&norm)) norm=true;
    /* if normalize */
    if (norm) {
	if (!sf_getfloat("clip",&clip)) clip = FLT_EPSILON;
	/* clip for fold normalization */
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);    

    if (interp) {
	count  = sf_floatalloc(nm);
	if (NULL != sf_getstring("fold")) {
	    /* output file for fold (optional) */ 
	    fold = sf_output("fold");
	    sf_putint(fold,"n1",nx);
	    sf_putint(fold,"n2",ny);
	    sf_putint(fold,"n3",n4);
	    sf_putfloat(fold,"o1",x0);
	    sf_putfloat(fold,"o2",y0);
	    sf_putfloat(fold,"d1",dx);
	    sf_putfloat(fold,"d2",dy);
	} else {
	    fold = NULL;
	}
    } else {
	count = NULL;
	fold = NULL;
    }

    for (i4=0; i4 < n4; i4++) {
	switch (interp) {
	    case 0:
		medbin_init (xy[i4], x0,y0,dx,dy,nx,ny, nd);
		sf_warning("Using median interpolation (i4= %d of %d);", i4+1,n4);
		break;
	    case 1:
		sf_int2_init (xy[i4], x0,y0,dx,dy,nx,ny, sf_bin_int, 1, nd);
		sf_warning("Using nearest-neighbor interpolation (i4= %d of %d);", i4+1,n4);
		break;
	    case 2:
		sf_int2_init (xy[i4], x0,y0,dx,dy,nx,ny, sf_lin_int, 2, nd);
		sf_warning("Using linear interpolation (i4= %d of %d);", i4+1,n4);
		break;
	    default:
		sf_error("Unsupported interp=%d",interp);
		break;
	}

	if (interp) { /* compute fold */
	    for (id=0; id<nd; id++) {
		dd[id]=1.;
	    }
	    sf_int2_lop (true, false,nm,nd,count,dd); 
	    if (NULL != fold) sf_floatwrite (count,nm,fold);
	
	    if (norm) {
		for (im=0; im<nm; im++) {
		    if (clip < count[im]) count[im]=1./fabsf(count[im]);
		    else                  count[im]=0.;
		}
	    }
	}

	for (it=0; it < nt; it++) { /* loop over time slices */
	    sf_floatread (dd,nd,in);
	    if (interp) {
		sf_int2_lop (true,false,nm,nd,mm,dd);
		if (norm) {
		    for (im=0; im<nm; im++) {
			mm[im] *= count[im];
		    }
		}
	    } else {
		medbin (dd,mm);
	    }
	    sf_floatwrite (mm,nm,out);
	}
    } /* i4 */
    sf_warning(".");


    exit(0);
}

/* 	$Id: Mbin.c 13815 2015-02-01 23:48:56Z sfomel $	 */
