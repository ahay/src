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

static float *t0;

static int fermat(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(int *)a];
    tb = t0[*(int *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}

int main (int argc, char *argv[])
{
    bool verb, up2, up3;
    unsigned char update;
    int n1,n2,n3, n12, n23, a2,b2,c2, a3,b3,c3, *visit, order, i,j,i1,i2,i3;
    float eps, **dat, **p, **q, *time, o1, d1, t2;
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

    if (NULL != sf_getstring("seed")) {
	seed = sf_input("seed");
	time = NULL;
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
    q = sf_floatalloc2(n1,n23);
    dat = sf_floatalloc2(n1,n23);

    sf_floatread(p[0],n12,dip);
    sf_floatread(q[0],n12,dip);

    t0 = sf_floatalloc(n23);
    visit = sf_intalloc(n23);
    sf_floatread(t0,n23,cost);

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init(n1,n2, eps*eps, order, 1, false);

    /* sort from small to large traveltime */
    for (i = 0; i < n23; i++) {
	visit[i] = i;
    }
    qsort(visit, n23, sizeof(int), fermat);

    for (i = 0; i < n23; i++) {
	j = visit[i];
	t2 = t0[j];

	i2 = j%n2;
	i3 = j/n2;

	update = 0;

	a2 = j-1;
	b2 = j+1;
	up2 = i2 && (i2 == n2-1 || 1 != fermat(&a2,&b2));
	c2 = up2? a2:b2;
	if (t2 < t0[c2]) update |= 1;

	a3 = j-n2;
	b3 = j+n2;
	up3 = i3 && (i3 == n3-1 || 1 != fermat(&a3,&b3));
	c3 = up3? a3:b3;
	if (t2 < t0[c3]) update |= 2;
	
	switch(update) {
	    case 0:
		if (NULL != seed) {
		    sf_floatread(dat[j],n1,seed);
		} else {
		    for (i1=0; i1 < n1; i1++) {
			dat[j][i1] = time[i1];
		    }
		}
		break;
	    case 1:
		predict1_step(up2,dat[c2],up2? p[a2]:p[j],dat[j]);
		break;
	    case 2:
		predict1_step(up3,dat[c3],up2? q[a3]:q[j],dat[j]);
		break;
	    case 3:
		predict2_step(up2,up3,dat[c2],dat[c3],
			      up2? p[a2]:p[j],
			      up3? q[a3]:q[j],dat[j]);
		break;
	}
    }
	
    sf_floatwrite(dat[0],n12,out);
    exit (0);
}
