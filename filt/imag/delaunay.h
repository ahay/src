/*
 * File: delaunay.h
 * ----------------
 * Interface for the incremental delaunay triangulation. 
 */
#ifndef _delaunay_h
#define _delaunay_h

#include <rsf.h>

#include "list_struct.h"

/* 
 * Function: DelaunayNew
 * ---------------------
 * Creates a big initial triangle, using the bounding box information.
 * Adds 3 nodes and 3 edges to the lists. 
 * The nodes are placed at 
 * {(3*xmin-xmax)/2, ymin}, {(3*xmax-xmin)/2}, 
 * {(xmax+xmin)/2, 2*ymax-ymin}
 */
void DelaunayNew (double xmin, double xmax, 
		  double ymin, double ymax, double zero);

/*
 * Function: FreeDelaunay
 * ----------------------
 * Frees the storage associated with the triangulation.
 * (Except nodes and edges.)
 */
void DelaunayFree (void);

/* 
 * Function: WriteDelaunay
 * -----------------------
 * Writes the nodes of the triangulation
 * to file. The output format is
 *
 * a_0 b_0 c_0
 * a_1 b_1 c_1
 * ...
 */
void DelaunayWrite (sf_file file);

/* 
 * Function: InsertNode
 * --------------------
 * Inserts node q to the triangulation.
 */
void InsertNode (Node q);

/* 
 * Function: InsertEdge;
 * --------------------
 * Inserts a boundary edge to the triangulation,
 * making the necessary adjustments.
 */ 
void InsertEdge (Edge ab);

/* 
 * Function: DelaunayRefine
 * ------------------------
 * Refine triangulation by inserting nr nodes.
 * Return the number left.
 */
int DelaunayRefine (int nr); 

/*
 * Function:Interpolate
 * --------------------
 * Find a value at a certain point
 * by linear triangular interpolation
 */
double Interpolate (Node q);

/*
 * Function:BoundingBox
 * --------------------
 * Initializes a bounding box
 */
void BoundingBox (double* box);

#endif
