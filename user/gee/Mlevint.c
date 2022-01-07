/* Leveler inverse interpolation in 1-D. */
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

#include "levint.h"

int main(int argc, char* argv[])
{
    int nd, m1, na, nr, niter;
    float *rr, *dd, *coord, o1, d1, eps;
    char *header;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* create model */
    if (!sf_getint ("nx",&m1)) sf_error("Need n1=");
    /* number of bins */
    sf_putint(out,"n1",m1);

    if (!sf_getfloat("x0",&o1)) sf_error("Need o1=");
    /* grid origin */
    sf_putfloat (out,"o1",o1);

    if (!sf_getfloat("dx",&d1)) sf_error("Need d1=");
    /* grid sampling */
    sf_putfloat (out,"d1",d1);

    if (!sf_getint("niter",&niter)) niter=1+m1*3/2; 
    /* number of conjugate-gradient iterations */
    niter *= 2;

    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* regularization parameter */

    /* create filter */
    if (!sf_getint("na",&na)) na=3;
    nr = m1 + na;

    rr = sf_floatalloc(nr);

    if (!sf_histint(in,"n1",&nd)) nd=1;
    coord = sf_floatalloc(nd);
    dd = sf_floatalloc(nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");
    sf_floatread (coord,nd,head);
    sf_fileclose (head);

    sf_floatread (dd,nd,in);

    levint1 (niter, m1, nr, nd, coord, dd, o1, d1, rr, eps);
    
    sf_floatwrite (rr,m1,out);


    exit(0);
}

/* 	$Id$	 */
