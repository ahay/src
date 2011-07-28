/* Incremental Delaunay triangulation. */
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
#include <string.h>
#include <math.h>

#include <rsf.h>
/*^*/

#include "delaunay.h"

#include "list_struct.h"
#include "_basic_struct.h"
/*^*/

#include "predicates.h"

/* private variables */
static Triangle BigTriangle;
static double* BBox = NULL;

/* private functions */
static double TestEdge (Edge ab, Node q);
static int GoodEdge (Edge ab);
static void InsertANode (Node q, Triangle abc);
static void InsertEdgeNode (Node q, Edge ab);
static void InsertInternalNode (Node q, Triangle abc);
static void GetChild (Triangle abc, int n);
static void Validate (Triangle abc);
static int ThisEdge (Triangle abc, Edge ab);
static Node Near (Edge ab, Triangle abc);
static Triangle LocateNode (Triangle abc, Node q); 
static void FreeTriangle (Triangle abc);
static void WriteTriangle (sf_file file, Triangle abc);
static Edge IsEdge (Triangle abc, Edge ab);
static Edge IsNode (Triangle abc, Edge ab);
static Node CrossNode (Node* one, Node* two);
static void Bounding (Node q);
static void Rivara (int* nr, Triangle abc);
static Edge LongestEdge (Triangle abc);
static double DotProd (Edge one, Edge two);
static void Mark (Triangle abc);

/* !!!! debugging !!!! 
static void InTriangle (Triangle abc) {
    WriteTriangle (stdout, abc);
    WriteNode (abc->edge[0]->ends[0]);
    WriteNode (abc->edge[0]->ends[1]);
    WriteNode (abc->edge[1]->ends[0]);
    WriteEdge (abc->edge[1]);
    WriteEdge (abc->edge[2]);
}
*/

void DelaunayNew (double xmin, double xmax, 
		  double ymin, double ymax, double zero)
/*< Creates a big initial triangle, using the bounding box information.
 * Adds 3 nodes  and 3 edges. 
 * The nodes are placed at 
 * {(3*xmin-xmax)/2, ymin}, {(3*xmax-xmin)/2}, 
 * {(xmax+xmin)/2, 2*ymax-ymin}
 >*/
{
    Edge ab;
    Node a[3];
    int i;

    a[0] = AppendNode ((3*xmin-xmax)/2, ymin,      zero, EMPTY);
    a[1] = AppendNode ((3*xmax-xmin)/2, ymin,      zero, EMPTY);
    a[2] = AppendNode ((xmax+xmin)/2, 2*ymax-ymin, zero, EMPTY);

    BigTriangle = (Triangle) malloc (sizeof (*BigTriangle));
    BigTriangle->child = NULL;
    BigTriangle->clone = BOUNDARY;

    BigTriangle->edge[0] = AppendEdge (a[0], a[1], EMPTY);
    BigTriangle->edge[0]->face[0] = BigTriangle;
    BigTriangle->edge[0]->face[1] = NULL;
    for (i = 1; i < 3; i++) {
	ab = BigTriangle->edge[i] = AppendEdge (a[2], a[2-i], EMPTY);
	ab->face[2-i] = BigTriangle;
	ab->face[i-1] = NULL;
    }
#ifndef FAST
    exactinit ();
#endif
}

static double TestEdge (Edge ab, Node q)
/* The output > 0 if q is on the left of ab */
{
    double test;

#ifdef FAST
    test = orient2dfast (ab->ends[0]->x,ab->ends[1]->x,q->x);
#else
    test = orient2d     (ab->ends[0]->x,ab->ends[1]->x,q->x);
#endif
    return test;
}

static int GoodEdge (Edge ab)
/* test if an edge is good */
{
    double test;
    Node near[2];

    if (ab->face[0] == NULL || ab->face[1] == NULL) return 1;

    near[0] = Near (ab, ab->face[0]);
    near[1] = Near (ab, ab->face[1]);


#ifdef FAST
    test = incirclefast (ab->ends[0]->x,near[0]->x,
			 ab->ends[1]->x,near[1]->x);
#else
    test = incircle     (ab->ends[0]->x,near[0]->x,
			 ab->ends[1]->x,near[1]->x);
#endif

    return (ab->type == BOUNDARY)? (test >= 0.) : (test > 0.);
}

static Triangle LocateNode (Triangle abc, Node q) 
/* Starting from the bounding triangle abc, outputs
 * the triangle, containing the node q. 
 */
{
    int i, k;

    if (abc->child == NULL) {
	if (BBox != NULL) {
	    Bounding (abc->edge[0]->ends[0]);
	    Bounding (abc->edge[0]->ends[1]);
	    Bounding (abc->edge[1]->ends[0]);
	}
	return abc;
    }
    if (abc->child[2] == NULL) { 
	k = (TestEdge (abc->child[0]->edge[2],q) < 0);
    } else {    
	for (i = 0, k = -1; i < 3; i++) {
	    if (TestEdge (abc->child[i]->edge[2],q) >= 0) {
		if (k >= 0) {
		    k = (5 - k - i)%3;
		    break;
		} else {
		    k = i;
		}
	    }
	}
    }
    return LocateNode (abc->child[k], q);
}

void DelaunayFree (void)
/*<  * Frees the storage associated with the triangulation.
 * (Except nodes and edges.)
 >*/
{
    FreeTriangle (BigTriangle);
}

static void FreeTriangle (Triangle abc)
{
    int i;

    if (abc->child != NULL && abc->clone != EMPTY) {
	for (i=0; abc->child[i] != NULL; i++) {
	    FreeTriangle (abc->child[i]);
	}
	free (abc->child);
    }
    free (abc);
}

void DelaunayWrite (sf_file file)
/*< Writes the nodes of the triangulation
 * to file. The output format is
 *
 * a_0 b_0 c_0
 * a_1 b_1 c_1
 >*/
{
    WriteTriangle (file, BigTriangle);
}
 
static void WriteTriangle (sf_file file, Triangle abc)
/* Writes the nodes of triangle abc and all its descendants
 * to file. 
 */
{
    int i, ends[3];

    if (abc->child != NULL) {
	if (abc->clone != EMPTY) {
	    for (i = 0; abc->child[i] != NULL; i++) {
		WriteTriangle (file, abc->child[i]);
	    }
	}
    } else {
	ends[0] = NodeNumber (abc->edge[0]->ends[0]);
	ends[1] = NodeNumber (abc->edge[0]->ends[1]);
	ends[2] = NodeNumber (abc->edge[1]->ends[0]);
	sf_intwrite(ends,3,file);
    }
}

static void InsertANode (Node q, Triangle abc)
{
    int i;
    Edge ab;

    for (i = 0, ab = NULL; i< 3; i++) { 
	if (TestEdge (abc->edge[i],q) == 0.) {
	    ab = abc->edge[i];
	    break;
	}
    }
       
    if (ab != NULL) {
	InsertEdgeNode (q, ab);
    } else {
	InsertInternalNode (q, abc);
    }
}

void InsertNode (Node q)
/*< Inserts node q to the triangulation. >*/
{
    if (BBox != NULL) {
	BBox[0] = BBox[1] = q->x[0];
	BBox[2] = BBox[3] = q->x[1];
    }
    InsertANode (q, LocateNode (BigTriangle, q));
}

static void InsertInternalNode (Node q, Triangle abc)
{
    int i, j;
    Triangle child;
    Edge ab, qa, qb;

    GetChild (abc, 3);
    qa = AppendEdge (q, abc->edge[1]->ends[0], ADDED);
    ab = abc->child[0]->edge[0] = abc->edge[0];
    j = (ab->face[0] != abc);
    ab->face[j] = abc->child[0];
    for (i = 1; i < 3; i++) {
	child = abc->child[i];
	child->edge[0]   = ab = abc->edge[i];
	child->edge[i]   = qa;
	qa->face[2-i] = child;
	child->edge[3-i] = qb = AppendEdge (q, ab->ends[1], ADDED);
	qb->face[i-1] = ab->face[2-i] = child;
	qb->face[2-i] = abc->child[0];
	abc->child[0]->edge[i] = qb;
    }
    for (i = 0; i < 3; i++) {
	Validate (abc->child[i]);
    }
}

static void InsertEdgeNode (Node q, Edge ab)
/* insert node to edge */
{
    int i, j, k, l;
    Triangle abc, *child;
    Edge qa;

    for (i = 0; i< 2; i++) {
	abc = ab->face[i];
	if (abc == NULL) fprintf (stderr,"Trouble\n");
	GetChild (abc, 2);
	child = abc->child;
	k = ThisEdge (abc, ab);
	for (l = 0; l < 2; l++) {
	    qa= child[l]->edge[0] = abc->edge[(k+2-l)%3];   
	    j = (qa->face[0] != abc);
	    qa->face[j] = child[l];
	}
	qa= child[0]->edge[2] = child[1]->edge[1] =
	    AppendEdge (q, Near (ab, abc), ADDED);
	qa->face[0] = child[0];
	qa->face[1] = child[1];
	qa = child[0]->edge[1] = 
	    AppendEdge (q, ab->ends[i], ADDED);
	qa->face[1] = child[0];
    }
    for (i = 0; i< 2; i++) {
	ab->face[i]->child[0]->edge[1]->face[0] = ab->face[!i]->child[1];
	ab->face[i]->child[1]->edge[2] = ab->face[!i]->child[0]->edge[1];
    }
    ab->type = EMPTY;
    for (i = 0; i< 2; i++) {
	abc = ab->face[i];
	Validate (abc->child[1]);
	Validate (abc->child[0]);
    }
}

static Node CrossNode (Node* one, Node* two)
{
    int i;
    double da[3], db[3], step;
    Node node;

    node = NULL;
    for (i=0; i<3; i++) {
	da[i] = one[1]->x[i] - one[0]->x[i];
	db[i] = two[1]->x[i] - two[0]->x[i];
    }
    step = da[0]*db[1]-da[1]*db[0];
    if (fabs (step) > 0.) { 
	step = (db[1]*(two[0]->x[0] - one[0]->x[0]) + 
		db[0]*(one[0]->x[1] - two[0]->x[1]))/
	    step;
	if (step > 0. && step < 1.) 
	    node = AppendNode (one[0]->x[0] + da[0]*step,
			       one[0]->x[1] + da[1]*step,
			       one[0]->x[2] + da[2]*step,
			       ADDED);
    }
    return node;
}

static void Validate (Triangle abc)
/* validate triangle */
{
    int i, j, k, l;
    Edge ab, qa, ac, ad;
    Node q, near[2];
    Triangle abd;
    Triangle* child;

    ab = abc->edge[0];
    if (GoodEdge (ab)) return;

    j = (ab->face[0] != abc);
    abd = ab->face[!j];
    near[0] = Near (ab, abc);
    near[1] = Near (ab, abd);
    if (BBox != NULL) Bounding (near[1]);

    if (ab->type == EMPTY && ab->ends[1]->type != EMPTY) {
	if ((q = CrossNode (ab->ends, near)) == NULL) {
	    fprintf (stderr,"Bad\n");
	} else {
	    InsertEdgeNode (q, ab);
	}
	return;
    }
    ab->type = EMPTY;
    abd->clone = EMPTY;
    GetChild (abc, 2);
    child = abd->child = abc->child;
    qa= AppendEdge (near[0], near[1], ADDED);
    k = ThisEdge (abd,ab);
    for (i=0; i<2; i++) {
	qa->face[i] = child[i];
	ac = child[i]->edge[0]   = abd->edge[(k+2-i)%3];
	ad = child[i]->edge[i+1] = abc->edge[i+1];
	l = (ac->face[0] != abd);
	ad->face[!i] = ac->face[l] = child[i];
	child[i]->edge[2-i] = qa;
    }
    Validate (child[1]);
    Validate (child[0]);
}

static int ThisEdge (Triangle abc, Edge ab)
{
    int i;
    for (i = 0; i < 3; i++) {
	if (ab == abc->edge[i]) return i;
    }
    return -1;
}

static void GetChild (Triangle abc, int n)
{
    int i;

    abc->child = (Triangle*) calloc (n+1,sizeof(Triangle));
    for (i=0; i < n; i++) {
	abc->child[i] = (Triangle) malloc (sizeof (*abc));
	abc->child[i]->child = NULL;
	abc->child[i]->clone = BOUNDARY;
    }  
    abc->child[n] = NULL;
}

static Node Near (Edge ab, Triangle abc) 
{
    int k;
    k = 3 - ThisEdge (abc,ab);
    return (k == 3)? abc->edge[1]->ends[0] : abc->edge[k]->ends[1];
}

void InsertEdge (Edge ab)
/*< Inserts a boundary edge to the triangulation,
 * making the necessary adjustments. >*/
{
    Triangle abc;
    Edge old;
    Node midpoint, newpoint;
    int i;

    midpoint = (Node) malloc (sizeof (*midpoint));
    midpoint->x = (double *) calloc (DIMENSION, sizeof(double));
    for (i=0; i<3; i++) {
	midpoint->x[i] = (ab->ends[0]->x[i]+
			  ab->ends[1]->x[i])/2;
    }
    midpoint->type = EMPTY;
  
    abc = LocateNode (BigTriangle, midpoint);
    free (midpoint->x);
    free (midpoint);

    if ((old = IsEdge (abc, ab)) != NULL) {
	old->type = EMPTY;
	return;
    }
    if ((old = IsNode (abc, ab)) != NULL) {
	newpoint = CrossNode (old->ends,ab->ends);
	if (newpoint == NULL) {
	    fprintf (stderr,"Too bad\n");
	    return;
	}
	if (old->face[0] == NULL || 
	    old->face[1] == NULL) {
	    fprintf (stderr,"Oops\n");
	    return;
	}
    } else {
	newpoint = NULL;
	for (i=0; i<3; i++) {
	    old = abc->edge[i];
	    if (old->face[0] == NULL || 
		old->face[1] == NULL) continue;
	    if ((newpoint = CrossNode (old->ends,ab->ends)) != NULL) break;
	}
	if (newpoint == NULL) {
	    fprintf (stderr,"Disaster\n");
	    return;
	}
    }
    InsertEdgeNode (newpoint,old);    
    InsertEdge (AppendEdge (newpoint, ab->ends[0],EMPTY));
    InsertEdge (AppendEdge (newpoint, ab->ends[1],EMPTY));
}

static Edge IsEdge (Triangle abc, Edge ab)
{
    int i;
    Edge test;

    if (abc->edge[1]->ends[0] == ab->ends[0]) {
	for (i=1; i < 3; i++) {
	    test = abc->edge[i];
	    if (test->ends[1] == ab->ends[1]) return test;
	}
    } else {
	test = abc->edge[0];
	for (i=0; i < 2; i++) {
	    if (test->ends[ i] == ab->ends[0] && 
		test->ends[!i] == ab->ends[1]) return test;
	}
    }
    return NULL;
}

static Edge IsNode (Triangle abc, Edge ab)
{
    int i;

    for (i=0; i<2; i++) {
	if (abc->edge[1]->ends[0] == ab->ends[i]) 
	    return abc->edge[0];
	if (abc->edge[1]->ends[1] == ab->ends[i])
	    return abc->edge[2];
	if (abc->edge[2]->ends[1] == ab->ends[i])
	    return abc->edge[1];
    }
    return NULL;
}

void BoundingBox (double* box)
/*< Initializes a bounding box >*/
{
    BBox = box;
}

static void Bounding (Node q)
{
    double x;

    if ((x = q->x[0]) < BBox[0]) {
	BBox[0] = x;
    } else if (x > BBox[1]) {
	BBox[1] = x;
    }

    if ((x = q->x[1]) < BBox[2]) {
	BBox[2] = x;
    } else if (x > BBox[3]) {
	BBox[3] = x;
    }
}

double Interpolate (Node q)
/*< Find a value at a certain point
 * by linear triangular interpolation >*/
{
    double *a, *b, *c, abc, xbc, axc, abx, z;
    Triangle t;

    t = LocateNode (BigTriangle, q);

    a = t->edge[0]->ends[0]->x;
    b = t->edge[0]->ends[1]->x;
    c = t->edge[1]->ends[0]->x;

    abc = orient2d (a,b,c);
    xbc = orient2d (q->x,b,c);
    axc = orient2d (a,q->x,c);
    abx = orient2d (a,b,q->x);

    z = (a[2]*xbc + b[2]*axc + c[2]*abx)/abc;
    return z;
}

int DelaunayRefine (int nr) 
/*< Refine triangulation by inserting nr nodes
 * Return the number left >*/
{
    Mark (BigTriangle);
    Rivara (&nr, BigTriangle);
    return nr;
}

static void Mark (Triangle abc)
{
    int i;
    Edge ab, ac, bc;

    if (abc->child != NULL) {
	if (abc->clone != EMPTY) {
	    for (i = 0; abc->child[i] != NULL; i++) {
		Mark (abc->child[i]);
	    }
	}
    } else if (abc->clone != ADDED) {
	ab = abc->edge[1]; if (ab->type == EMPTY) return;
	ac = abc->edge[2]; if (ac->type == EMPTY) return;
	bc = abc->edge[0]; if (bc->type == EMPTY) return;
	if ((DotProd (ab,ac) < 0.75*EdgeLength(ab)*EdgeLength(ac)) &&
	    (DotProd (ab,bc) < 0.75*EdgeLength(ab)*EdgeLength(bc)) &&
	    (DotProd (ac,bc) < 0.75*EdgeLength(ac)*EdgeLength(bc))) return;
	abc->clone = ADDED;
    }
    return;
}

static void Rivara (int* nr, Triangle abc)
/* Rivara refinement algorithm */
{
    int i;
    Edge ab;
    double x[3];

    if (abc->child != NULL) {
	if (abc->clone != EMPTY) {
	    for (i = 0; abc->child[i] != NULL; i++) {
		Rivara (nr, abc->child[i]);
	    }
	}
    } else if (abc->clone == ADDED) {
	while ((*nr > 0) && (abc->child == NULL)) {
	    ab = LongestEdge (abc);
	    for (i=0; i<3; i++) {
		x[i] = 0.5*(ab->ends[0]->x[i]+
			    ab->ends[1]->x[i]);
	    }
	    InsertEdgeNode (AppendNode(x[0],x[1],x[2],ADDED), ab);
	    (*nr)--;
	}
    }
    return;
}

static double DotProd (Edge one, Edge two)
/* dot product */
{
    int i;
    double dotprod;
  
    dotprod=0.;
    for (i=0; i <2; i++) {
	dotprod += 
	    (one->ends[1]->x[i] - one->ends[0]->x[i]) *
	    (two->ends[1]->x[i] - two->ends[0]->x[i]);
    }
    return (dotprod*dotprod);
}

static Edge LongestEdge (Triangle abc)
{
    Triangle abd;
    Edge ab, edge;
    int i, j, k;
    double length, lab;

    edge = NULL;
    for (length=0.,i=0;i<3;i++) {
	ab = abc->edge[i];
	if ((ab != edge) && 
	    (ab->type != EMPTY) && 
	    (length < (lab = EdgeLength(ab)))) {
	    edge = ab;
	    length = lab;
	}
    }
  
    while (1) {
	j = (edge->face[0] != abc);
	abd = edge->face[!j];
	for (k=i=0;i<3;i++) {
	    ab = abd->edge[i];
	    if ((ab != edge) && 
		(ab->type != EMPTY) && 
		(length < (lab = EdgeLength(ab)))) {
		edge = ab;
		length = lab;
		k++;
	    }
	}   
	if (k == 0) return edge;
	abc = abd;
    }
}

