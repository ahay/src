/* Decimate a regular grid to triplets for triangulation. */
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

#include "list_struct.h"
#include "delaunay.h"
#include "heap.h"

int main(int argc, char* argv[])
{
    int n1, n2, n12, i1, i2, i, it, nt, ne, b1, e1, b2, e2;
    float o1, o2, d1, d2, zero, **data=NULL, **xyz=NULL, *e=NULL;
    double BBox[4], error;
    hpoint *heap=NULL, *h=NULL;
    sf_file in=NULL, out=NULL, edge=NULL;
    Node q;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* number of triplets */

    if (!sf_getfloat("zero",&zero)) zero = 0.;
    /* level surface */

    sf_putint(out,"n1",3);
    sf_putint(out,"n2",nt);

    n12 = n1*n2;
    heap = (hpoint *) sf_alloc (n12,sizeof(hpoint));

    data = sf_floatalloc2(n1,n2);
    sf_floatread(data[0],n12,in);

    CreateNodeList (7+n12);
    CreateEdgeList ();

    heap_init (n12);

    for (i2 =0; i2 < n2; i2++) {
	for (i1 =0; i1 < n1; i1++) {
	    AppendNode ((double) (o1+i1*d1), 
			(double) (o2+i2*d2), 
			(double) data[i2][i1], BOUNDARY);
	}
    }

    DelaunayNew ((double) (o1-2*n1*d1), (double) (o1+2*n1*d1), 
		 (double) (o2-2*n2*d2), (double) (o2+2*n2*d2), (double) zero);
    BoundingBox (BBox);

    xyz = sf_floatalloc2(3,nt);
 
    InsertNode (AppendNode ((double) (o1-   d1), 
			    (double) (o2-   d2), 
			    (double) data[0][0],      BOUNDARY));
    InsertNode (AppendNode ((double) (o1+n1*d1), 
			    (double) (o2-   d2), 
			    (double) data[0][n1-1],   BOUNDARY));
    InsertNode (AppendNode ((double) (o1-   d1), 
			    (double) (o2+n2*d2), 
			    (double) data[n2-1][0], BOUNDARY));
    InsertNode (AppendNode ((double) (o1+n1*d1), 
			    (double) (o2+n2*d2), 
			    (double) data[n2-1][n1-1],  BOUNDARY));

    for (i=0; i <n12; i++) {
	i2 = i%n1;
	i1 = i-i2*n1;
	h = heap + i;
	error = (double) data[i2][i1] - Interpolate (GetNode(i));
	h->v = error*error;
	heap_insert (h);
    }

    for (it=0; it < nt; it++) {
	h = heap_extract();
	i = h - heap;
	q = GetNode(i);
	NodeOut (q,xyz[it],xyz[it]+1,xyz[it]+2);
	InsertNode (q);
	b1 = (BBox[0]-o1)/d1; e1 = (BBox[1]-o1)/d1;
	b2 = (BBox[2]-o2)/d2; e2 = (BBox[3]-o2)/d2;
	if (b1 < 0) b1 = 0; if (e1 > n1) e1 = n1;
	if (b2 < 0) b2 = 0; if (e2 > n2) e2 = n2;
	for (i2 = b2; i2 < e2; i2++) {
	    for (i1 = b1; i1 < e1; i1++) {
		i = i1+i2*n1;
		h = heap + i;
		error = (double) data[i2][i1] - Interpolate (GetNode(i));
		h->v = error*error;
		heap_update (h);
	    }
	}
    }

    sf_floatwrite(xyz[0],nt*3,out);

    if (NULL != sf_getstring("edgeout")) {
	e = sf_floatalloc (100*nt);
	ne = EdgeOut (e);

	edge = sf_output("edgeout");
	sf_putint(edge,"n1",2);
	sf_putint(edge,"n2",ne);
	sf_settype(edge,SF_COMPLEX);
	sf_fileflush(edge,in);

	sf_settype(edge,SF_FLOAT);
	sf_floatwrite(e,4*ne,edge);
    }

    exit (0);
}
