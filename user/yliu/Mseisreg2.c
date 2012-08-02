/* Data regularization in 2-D using seislet transform. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <rsfpwd.h>

#include "int1.h"

int main(int argc, char* argv[])
{
    int niter, nw, n1, i3, n3, nt, nd, nm, interp, order;
    float *mm, *dd, **pp, *offset, x0, dx, eps;
    char *header;
    char *type;
    bool verb;
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
    offset = sf_floatalloc(nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");
    sf_floatread (offset,nd,head);
    sf_fileclose (head);
    
    /* regular data */
    if (!sf_histint(dip,"n1",&n1) || n1 != nt) sf_error("Need n1=%d in dip",nt);
    if (!sf_histint(dip,"n2",&nm)) sf_error("Need n2= in dip");
    if (SF_FLOAT != sf_gettype(dip)) sf_error("Need float data in dip");

    if (!sf_histfloat(dip,"o2",&x0)) sf_error("Need o2= in dip");
    if (!sf_histfloat(dip,"d2",&dx)) sf_error("Need d2= in dip");

    sf_putint(out,"n2",nm);
    sf_putfloat(out,"o2",x0);
    sf_putfloat(out,"d2",dx);

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    if (!sf_getint("niter",&niter)) niter=50;
    /* number of iterations */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    pp = sf_floatalloc2(nt,nm);
    mm = sf_floatalloc(nt*nm);
    dd = sf_floatalloc(nt*nd);

    int1_init (offset, x0,dx,nm, nt, sf_spline_int, interp, nd);
    seislet_init(nt,nm,true,true,eps,order,type[0]);
    seislet_set(pp);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	/* read irregular data */
	sf_floatread(dd,nt*nd,in);

	/* read first dip */
	sf_floatread(pp[0],nt*nm,dip);

	sf_solver_reg(int1_lop, sf_cgstep, seislet_destruct, nt*nm, nt*nm, 
		      nt*nd, mm, dd, niter, eps, "verb", verb, "end");
 
	sf_cgstep_close();
	sf_floatwrite (mm,nt*nm,out);
    }


    exit(0);
}

/* 	$Id$	 */
