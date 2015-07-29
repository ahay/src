/* Manipulating structures for triangulation. */
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

#include <string.h>
#include <stddef.h>

#include <rsf.h>

#include "list_struct.h"
#include "_basic_struct.h"

#ifndef _list_struct_h

#define DIMENSION 3
/* number of dimensions */
/*^*/

enum elemType {EMPTY = -1,BOUNDARY,ADDED};
/* Flag to distinguish elements */
/*^*/

typedef struct CEdge *Edge;
typedef struct CNode *Node;
/* abstract data types */
/*^*/

#endif

/* private variables */
static Node NodeList, LastNode;
static Edge EdgeList, LastEdge;
static int NumbNode;

/* private functions */
static void FreeEdge (Edge edge);

void CreateNodeList (int n) 
/*< Allocates an empty list >*/
{
    LastNode = NodeList = (Node) sf_alloc (n,sizeof (CNode));
    NumbNode = n;
}   

void FreeNodeList (void) 
/*< Frees all the storage associated with the list. >*/
{
    Node  node;
    for (node = NodeList; node != LastNode; node++) {
	free (node->x);
    }
    free (NodeList);
}

Node AppendNode (double x, double y, double z, enum elemType type) 
/*< Appends a node to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.) >*/
{  
    Node node;

    if (LastNode == NodeList + NumbNode) return NULL; 
    node = LastNode;

    node->x = (double *) sf_alloc (DIMENSION, sizeof(double));
    node->x[0] = x;
    node->x[1] = y;
    node->x[2] = z;
    node->type = type;

    LastNode++;
    return node;
} 

int NodeNumber (Node node) 
/*< Translates a pointer to a node into the node number in the list.
 * If node == NULL, returns the total number of nodes. >*/
{
    int n;

    if (node == NULL) node = LastNode;
    n = node - NodeList;
    return n;
}

Node GetNode (int n)
/*< Returns n-th node. >*/
{
    return (NodeList+n);
}

void ReadEdgeList (sf_file edgefile)
/*< Reads an edge list from the file edgefile. >*/
{
    int i, lr[2], edges, two;

    if (SF_INT != sf_gettype(edgefile)) 
	sf_error("%s: Need int input for edges",__FILE__);
    if (!sf_histint(edgefile,"n1",&two) || 2 != two)
	sf_error("%s: Need n1=2 for edges",__FILE__);
    if (!sf_histint(edgefile,"n2",&edges)) edges=1; 
	
    for(i = 0; i < edges; i++) {
	sf_intread(lr,2,edgefile);
	if (lr[1] > lr[0]) { 
	    AppendEdge (NodeList+lr[1], NodeList+lr[0], BOUNDARY);
	} else {
	    AppendEdge (NodeList+lr[0], NodeList+lr[1], BOUNDARY);
	}
    }
}

void WriteList (sf_file nodefile, sf_file edgefile)
/*< Writes a  node list to the file nodefile.
 * Writes an edge list to the file edgefile. >*/
{
    Node node;
    Edge edge;
    float x[3];
    int ends[2];

    sf_putint(nodefile,"n1",3);
    sf_putint(nodefile,"n2",(int) (LastNode - NodeList));
    sf_settype(nodefile,SF_FLOAT);

     for (node = NodeList; node != LastNode; node++) {
	 x[0] = node->x[0];
	 x[1] = node->x[1];
	 x[2] = node->x[2];
	 sf_floatwrite(x,3,nodefile);
     }

     sf_putint(edgefile,"n1",2);
     sf_putint(edgefile,"n2",EdgeNumber ());
     sf_settype(edgefile,SF_INT);

    for (edge = EdgeList; edge->next != NULL; edge = edge->next) {
	if (edge->type != EMPTY) {
	    ends[0] = NodeNumber (edge->ends[0]);
	    ends[1] = NodeNumber (edge->ends[1]);
	    sf_intwrite(ends,2,edgefile);
	}
    }
}

void NodeExec (void (* exec) (Node))
/*< Recursive operations on nodes in the list. >*/
{
    Node node;

    for (node = NodeList; node != LastNode; node++) {
	if (node->type != EMPTY) {
	    exec (node);
	}
    } 
}

void NodeOut (Node q, float *x, float *y, float *z)
/*< Copy a Nodes to 3 numbers. >*/
{
    *x = q->x[0];
    *y = q->x[1];
    *z = q->x[2];
}

int EdgeOut (float *e)
/*< Copy Edges to a float array. Returns the number of edges. >*/
{
    int i, j;
    Edge edge;

    i = j = 0;
    for (edge = EdgeList; edge->next != NULL; edge = edge->next) {
	if (edge->type == BOUNDARY) {
	    e[i] = edge->ends[0]->x[0]; i++;
	    e[i] = edge->ends[0]->x[1]; i++; 
	    e[i] = edge->ends[1]->x[0]; i++;
	    e[i] = edge->ends[1]->x[1]; i++;
	    j++;
	}
    }
    for (edge = EdgeList; edge->next != NULL; edge = edge->next) {
	if (edge->type == ADDED) {
	    e[i] = edge->ends[0]->x[0]; i++;
	    e[i] = edge->ends[0]->x[1]; i++; 
	    e[i] = edge->ends[1]->x[0]; i++;
	    e[i] = edge->ends[1]->x[1]; i++;
	    j++;
	}
    }
    return j;
}

void WriteNode (Node list)
/*< Debugging >*/
{
    printf ("%lf %lf %lf %i\n", 
	    list->x[0], list->x[1], list->x[2], list->type);
}

void WriteEdge (Edge list)
/*< Debugging >*/
{
    printf ("%u %u %i\n", 
	    NodeNumber (list->ends[0]), 
	    NodeNumber (list->ends[1]), list->type);
    WriteNode (list->ends[0]);
    WriteNode (list->ends[1]);
}

void CreateEdgeList (void)
/*< Allocates an empty list >*/
{
    int i;

    LastEdge = EdgeList = (Edge) malloc (sizeof (CEdge));
    for (i=0 ; i < 2; i++) {
	EdgeList->ends[i] = NULL;
	EdgeList->face[i] = NULL;
    }
    EdgeList->next = NULL;
    EdgeList->type = BOUNDARY;
}   

void ClearEdgeList (void)
/*< Frees added edges in the list. >*/
{
    Edge tmp;
  
    while (EdgeList->type != BOUNDARY) {
	tmp = EdgeList;
	EdgeList = EdgeList->next;
	free (tmp);
    }
    if (EdgeList->next == NULL) {
	FreeEdgeList ();
	CreateEdgeList ();
	return;
    }
    LastEdge = EdgeList; 
    while (LastEdge->next->next != NULL) { 
	if (LastEdge->next->type == BOUNDARY) {
	    LastEdge = LastEdge->next;
	} else {
	    tmp = LastEdge->next;
	    LastEdge->next = tmp->next;
	    free (tmp);
	}
    }
    LastEdge = LastEdge->next;
}

void FreeEdgeList (void)
/*< Frees all the storage associated with the list. >*/
{
    FreeEdge (EdgeList);
}

static void FreeEdge (Edge edge)
/* Frees all the storage associated with 
 * the edge and its descendants recursively. */
{
    if (edge != NULL) {
	FreeEdge (edge->next);
	free (edge);
    }
}

Edge AppendEdge (Node left, Node right, enum elemType type) 
/*< Appends an edge to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.) >*/
{
    int i;
    Edge tmp;

    tmp = LastEdge;
    tmp->ends[0] = left;
    tmp->ends[1] = right;
    tmp->type = (right->type == EMPTY)? EMPTY : type;

    LastEdge = tmp->next = (Edge) malloc (sizeof (CEdge));
    for (i = 0; i < 2; i++) {
	LastEdge->ends[i] = NULL;
	LastEdge->face[i] = NULL;
    }
    LastEdge->type = BOUNDARY;
    LastEdge->next = NULL;

    return tmp;
} 

int EdgeNumber (void)
/*< Returns the number of edges in the list. >*/
{
    int n;
    Edge list;

    for (n = 0, list = EdgeList; list->next != NULL; list =list->next) { 
	if (list->type != EMPTY) n++;
    }
    return n; 
}

void EdgeEnumerate (void)
/*< Give numbers to edges. >*/
{
    int n;
    Edge list;

    for (n = 0, list = EdgeList; list->next != NULL; list =list->next) { 
	if (list->type != EMPTY) {
	    list->num = n;
	    list->type = BOUNDARY;
	    n++;
	}
    }
}

void EdgeExec (void (* exec) (Edge))
/*< Recursive operations on edges in the list. >*/
{
    Edge edge;

    for (edge = EdgeList; edge->next != NULL; edge = edge->next) {
	if (edge->type == BOUNDARY) {
	    exec (edge);
	}
    } 
}

double EdgeLength (Edge ab)
/*< Returns the length of an edge. >*/
{
    double dx, dy;

    dx = ab->ends[1]->x[0] - ab->ends[0]->x[0];
    dy = ab->ends[1]->x[1] - ab->ends[0]->x[1];
    return (dx*dx + dy*dy);
}

void MoveNode (Node q, double x, double y)
/*< New position for the node. >*/
{
    q->x[0] = x;
    q->x[1] = y;
}

void NodeValues(int s, int n, const float *value) {
/*< go through the list of nodes starting from s and reset values >*/
    int i;
    Node node;

    if (s < 0 || NodeList+s+n > LastNode) sf_error("%s: out of range",__FILE__);

    for (i=0; i < n; i++) {
	node = NodeList+s+i;
	node->x[2]=value[i];
    }
}
