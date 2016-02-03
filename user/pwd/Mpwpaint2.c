/* 3-D painting by plane-wave construction. */
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
#include <math.h>

#include <rsf.h>

#include "predict.h"
#include "update.h"

int main (int argc, char *argv[])
{
    bool verb, up2, up3, data;
    unsigned char update;
    int n1,n2,n3, n12, n23, j2, j3=0, order, i,j,i1;
    float eps, **dat, **p, **q, *p2, *p3=NULL, *time, *t0, o1, d1;
    sf_file dip, out, seed, cost;

    sf_init(argc,argv);
    dip = sf_input("in");
    out = sf_output("out");
    cost = sf_input("cost");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(dip,"n3",&n3)) sf_error("No n3= in input");
    n23 = n2*n3;
    n12 = n1*n23;
    sf_putint(out,"n4",1);

    if (!sf_getbool("data",&data)) data=false;
    /* spray input data */
    
    if (NULL != sf_getstring("seed")) {
	seed = sf_input("seed");
	time = NULL;
    } else if (data) {
	sf_error("Need seed=");
    } else {
	seed = NULL;
	time = sf_floatalloc(n1);
	if (!sf_histfloat(dip,"o1",&o1)) o1=0.; 
	if (!sf_histfloat(dip,"d1",&d1)) d1=1.; 
	for (i1=0; i1 < n1; i1++) {
	    time[i1] = o1+i1*d1;
	}
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    p = sf_floatalloc2(n1,n23);
    if (n3 > 1) {
	q = sf_floatalloc2(n1,n23);
    } else {
	q = NULL;
    }
    dat = sf_floatalloc2(n1,n23);

    sf_floatread(p[0],n12,dip);
    if (n3 > 1) sf_floatread(q[0],n12,dip);

    t0 = sf_floatalloc(n23);
    sf_floatread(t0,n23,cost);

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init(n1,n2, eps*eps, order, 1, true);

    update_init(n2,n3,t0);

    if (data) sf_floatread(dat[0],n12,seed);

    for (i = 0; i < n23; i++) {
	update = get_update(i,&up2,&up3,&j);
	
	j2 = up2? (j-1) :(j+1);
	if (NULL != q) j3 = up3? (j-n2):(j+n2);

	p2 = up2? p[j-1]: p[j];
	if (NULL != q) p3 = up3? q[j-n2]:q[j];

	switch(update) {
	    case 0:
		if (!data) {
		    if (NULL != seed) {
			sf_floatread(dat[j],n1,seed);
		    } else {
			for (i1=0; i1 < n1; i1++) {
			    dat[j][i1] = time[i1];
			}
		    }
		}
		break;
	    case 1:
		predict1_step(up2,dat[j2],p2,dat[j]);
		break;
	    case 2:
		if (NULL == q) sf_error("Wrong case");
		
		predict1_step(up3,dat[j3],p3,dat[j]);
		break;
	    case 3:
		if (NULL == q) sf_error("Wrong case");
		
		predict2_step(up2,up3,dat[j2],dat[j3],
			      p2,p3,dat[j]);
		break;
	}
    }
	
    sf_floatwrite(dat[0],n12,out);
    exit (0);
}
