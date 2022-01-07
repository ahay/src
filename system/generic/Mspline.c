/* 1-D cubic spline interpolation.

Specify either n1= o1= d1= or pattern=
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

#include <rsf.h>

#include "spline3.h"

static float **table=NULL, **table2=NULL;

static int compare_table (const void *a, const void *b)
{
    float fa = table[* ((int*) a)][0];
    float fb = table[* ((int*) b)][0];
 
    return (fa < fb)? -1: (fa > fb)? 1: 0;
}

int main(int argc, char* argv[])
{
    int nd, n1, i1, n2, i2, two, id, *index=NULL;
    float x, o1, d1, *table1=NULL, *trace=NULL, x0, dx, *fp=NULL;
    bool reginput, sort;
    sf_file in=NULL, out=NULL, pattern=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&two)) sf_error("Need n1= in input");
    if (two != 2) {
	reginput = true;
	nd = two;
	if (!sf_histfloat(in,"d1",&dx)) sf_error("Need d1= in input");
	if (!sf_histfloat(in,"o1",&x0)) sf_error("Need o1= in input");
	n2 = sf_leftsize(in,1);

	fp = NULL;
    } else {
	reginput = false;
	if (!sf_getbool("sort",&sort)) sort=false;
	/* if y, the coordinates need sorting */

	if (!sf_histint(in,"n2",&nd)) sf_error ("Need n2= in input");
	sf_putint(out,"n2",1);
	n2 = sf_leftsize(in,2);
	
	fp = sf_floatalloc(2);
	if (!sf_getfloats("fp",fp,2)) {
	    /* end-point derivatives */
	    free(fp);
	    fp = NULL;
	}
    }

    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }

    if (!sf_getint("n1",&n1) && 
	(NULL== pattern ||
	 !sf_histint(pattern,"n1",&n1))) sf_error("Need n1=");
    /* Output grid size */
    if (!sf_getfloat("d1",&d1) && 
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d1",&d1))) sf_error("Need d1=");
    /* Output sampling */
    if (!sf_getfloat("o1",&o1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o1",&o1))) sf_error("Need o1=");
    /* Output origin */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);

    trace = sf_floatalloc(n1);

    if (reginput) {
	spline3_init1(nd,x0,dx);
	table1 = sf_floatalloc(nd);
    } else {
	spline3_init(nd,fp);
	table = sf_floatalloc2(2,nd);
	if (sort) {
	  index = sf_intalloc(nd);
	  table2 = sf_floatalloc2(2,nd);
	}
    }

    for (i2=0; i2 < n2; i2++) {
	if (reginput) {
	    sf_floatread(table1,nd,in);
	    spline_coeffs1(table1);

	    for (i1=0; i1 < n1; i1++) {
		x = o1 + i1*d1;
		trace[i1] = spline_eval1(x);
	    }
	} else {
	    sf_floatread(table[0],2*nd,in);

	    if (sort) {
	      for (id=0; id < nd; id++) {
		index[id]=id;
	      }
	      qsort(index,nd,sizeof(int),compare_table);
	      for (id=0; id < nd; id++) {
		table2[id][0] = table[index[id]][0];
		table2[id][1] = table[index[id]][1];
	      }

	      spline_coeffs(table2);
	    } else {
	      spline_coeffs(table);
	    }

	    for (i1=0; i1 < n1; i1++) {
		x = o1 + i1*d1;
		trace[i1] = spline_eval(x);
	    }
	}

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
