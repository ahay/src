/* Data regularization in 2-D using plane-wave destruction. */
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
#include "twoplane2.h"
#include "allp2.h"
#include "int1.h"

int main(int argc, char* argv[])
{
    int niter, nw, n1, np, i3, n3, nt, nd, nm, interp;
    float *mm, *dd, **pp, **qq, *offset, x0, dx, eps;
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

    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    int1_init (offset, x0,dx,nm, nt, sf_spline_int, interp, nd);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_histint(dip,"n3",&np)) np=1;

    pp = sf_floatalloc2(nt,nm);

    if (np > 1) {
	qq = sf_floatalloc2(nt,nm);
    } else {
	qq = NULL;
    }

    mm = sf_floatalloc(nt*nm);
    dd = sf_floatalloc(nt*nd);

    if (NULL != qq) {
	twoplane2_init(nw, 1,1, nt,nm, drift, pp, qq);
    } else {
	allpass22_init(allpass2_init(nw, 1, nt,nm, drift, pp));
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */
 
    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	/* read irregular data */
	sf_floatread(dd,nt*nd,in);

	/* read first dip */
	sf_floatread(pp[0],nt*nm,dip);

	if (NULL != qq) {
	    /* read second dip */
	    sf_floatread(qq[0],nt*nm,dip);

	    sf_solver_reg(int1_lop, sf_cgstep, twoplane2_lop, nt*nm, nt*nm, nt*nd, mm, dd, niter, eps, "verb", verb, "end");
	} else {
	    sf_solver_reg(int1_lop, sf_cgstep, allpass21_lop, nt*nm, nt*nm, nt*nd, mm, dd, niter, eps, "verb", verb, "end");
	}
	sf_cgstep_close();
	
	/* write regular model */
	sf_floatwrite (mm,nt*nm,out);
    }
	


    exit(0);
}
