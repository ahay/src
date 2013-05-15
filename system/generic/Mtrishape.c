/* 2-D irregular data interpolation using triangulation and shaping regularization. */
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

#include <stdlib.h>

#include <rsf.h>

#include "list_struct.h"
#include "delaunay.h"

static int nd, n1, n2;
static float o1,d1, o2,d2;
static float *d;
static Node q;
static sf_triangle tr1=NULL, tr2=NULL;

static void shape(int n12, const float *inp, float *out, void *data)
/* I + S (BF - I) */
{
    int i1, i2, i;

    /* B F - I */

    sf_int2_lop (false,false,n12,nd,(float*) inp,d);
    NodeValues(3, nd, d);

    for (i=0; i < n12; i++) {
	i1 = i%n1;
	i2 = i/n1;

	MoveNode (q,  o1+i1*d1,  o2+i2*d2);
	out[i] = Interpolate (q) - inp[i];
    }

    /* S */

    if (NULL != tr1) {
	for (i2=0; i2 < n2; i2++) {
	    sf_smooth2 (tr1, 0, 1, false, false, out+i2*n1);
	}
    }
    if (NULL != tr2) {
	for (i1=0; i1 < n1; i1++) {
	    sf_smooth2 (tr2, i1, n1, false, false, out);
	}
    }   

    /* + I */

    for (i=0; i < n12; i++) {
	out[i] += inp[i];
    }
}
    
int main(int argc, char* argv[])
{
    float g1, g2, o3, g3;
    int n12, id, i2, i1, three, iter, niter, rect1, rect2, nw;
    float **xyz, **z, **m;
    float zero, xi, xmax, xmin, ymin, ymax, dx, dy, dz;
    sf_file in, out, pattern;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&three) || 3 != three) 
	sf_error("Need n1=3 in input");
    if (!sf_histint(in,"n2",&nd)) sf_error("No n2= in input");

    if (NULL != sf_getstring("pattern")) {
	/* pattern file for output dimensions */
	pattern = sf_input("pattern");
	
	if (!sf_histint(pattern,"n1",&n1)) sf_error("No n1= in pattern");
	if (!sf_histint(pattern,"n2",&n2)) sf_error("No n2= in pattern");
	if (!sf_histfloat(pattern,"d1",&d1)) d1=1.;
	if (!sf_histfloat(pattern,"d2",&d2)) d2=1.;
	if (!sf_histfloat(pattern,"o1",&o1)) o1=0.;
	if (!sf_histfloat(pattern,"o2",&o2)) o2=0.;
	
	sf_fileclose(pattern);
    } else {
	if (!sf_getint("n1",&n1)) sf_error("Need n1=");
	if (!sf_getint("n2",&n2)) sf_error("Need n2=");
	if (!sf_getfloat("d1",&d1)) d1=1.;
	if (!sf_getfloat("d2",&d2)) d2=1.;
	if (!sf_getfloat("o1",&o1)) o1=0.;
	if (!sf_getfloat("o2",&o2)) o2=0.;
    }

    n12 = n1*n2;

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);

    if (!sf_getfloat("zero",&zero)) zero = 0.;
    /* level surface */

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */

    if (niter > 0) {
	if (!sf_getint("rect1",&rect1)) rect1=1;
	if (!sf_getint("rect2",&rect2)) rect2=1;
	/* smoothing regularization */

	if (rect1 > 1) tr1 = sf_triangle_init (rect1,n1);
	if (rect2 > 1) tr2 = sf_triangle_init (rect2,n2);
    }

    if (!sf_getint("nw",&nw)) nw=2;
    /* interpolator size */
    
    xyz = sf_floatalloc2(3,nd);
    sf_floatread(xyz[0],nd*3,in);

    xmax = xmin = xyz[0][0]; 
    ymax = ymin = xyz[0][1]; 
    o3 = g3 = xyz[0][2];
    for (id =0; id < nd; id++) {
	if ((xi = xyz[id][0]) > xmax) {
	    xmax = xi;
	} else if (xi < xmin) {
	    xmin = xi;
	}
	if ((xi = xyz[id][1]) > ymax) {
	    ymax = xi;
	} else if (xi < ymin) {
	    ymin = xi;
	}
	if ((xi = xyz[id][2]) > g3) {
	    g3 = xi;
	} else if (xi < o3) {
	    o3 = xi;
	}    
    }
    dx = xmax - xmin; 
    dy = ymax - ymin;
    dz = g3-o3;

    xmin = xmin - dx; xmax = xmax + dx;
    ymin = ymin - dy; ymax = ymax + dy;
    o3 = o3 - dz;
    g3 = g3 + dz;

    sf_putfloat(out,"xmin",xmin);
    sf_putfloat(out,"xmax",xmax);
    sf_putfloat(out,"ymin",ymin);
    sf_putfloat(out,"ymax",ymax);

    g1 = o1+(n1-1)*d1;
    g2 = o2+(n2-1)*d2;

    if (o1 < xmin || g1 > xmax) sf_error("frame1 is too large\n");
    if (o2 < ymin || g2 > ymax) sf_error("frame2 is too large\n");

    CreateNodeList (4+nd);
    CreateEdgeList ();
    DelaunayNew ((double) xmin, (double) xmax, 
		 (double) ymin, (double) ymax, (double) zero);

    for (id =0; id < nd; id++) {
	InsertNode(AppendNode ((double) xyz[id][0], 
			       (double) xyz[id][1], 
			       (double) xyz[id][2], BOUNDARY));
    }

    z = sf_floatalloc2 (n1,n2);

    if (niter > 0) {
	m = sf_floatalloc2 (n1,n2);
	d = sf_floatalloc(nd);
    } else {
	m = z;
	d = NULL;
    }

    q = AppendNode (0.,0.,0.,EMPTY);
    for (i2 =0; i2 < n2; i2++) {
	for (i1 =0; i1 < n1; i1++) {
	    MoveNode (q,  o1+i1*d1,  o2+i2*d2);
	    z[i2][i1] = Interpolate (q);
	}
    }

    if (niter > 0) {
	sf_int2_init (xyz, o1, o2, d1, d2, n1, n2, sf_lg_int, nw, nd);

	for (i2 =0; i2 < n2; i2++) {
	    for (i1 =0; i1 < n1; i1++) {
		m[i2][i1] = z[i2][i1];
	    }
	}

	for (iter=0; iter < niter; iter++) {
	    /* m -> d */
	    sf_int2_lop (false,false,n12,nd,m[0],d);
	    NodeValues(3, nd, d);

	    for (i2 =0; i2 < n2; i2++) {
		for (i1 =0; i1 < n1; i1++) {
		    MoveNode (q,  o1+i1*d1,  o2+i2*d2);
		    m[i2][i1] += z[i2][i1] - Interpolate (q);
		}
	    }
	    
	    /* 2-D smoothing */
	    if (NULL != tr1) {
		for (i2=0; i2 < n2; i2++) {
		    sf_smooth2 (tr1, 0, 1, false, false, m[i2]);
		}
	    }
	    if (NULL != tr2) {
		for (i1=0; i1 < n1; i1++) {
		    sf_smooth2 (tr2, i1, n1, false, false, m[0]);
		}
	    }   
	}
    }
    
    sf_floatwrite (m[0], n12, out);

    exit(0);
}
