/*
 * File: list_struct.c
 * -------------------
 * Implementation for node and edge list manipulations
 */
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include <rsf.h>

#include "list_struct.h"
#include "_basic_struct.h"

/* private variables */
static Node NodeList, LastNode;
static Edge EdgeList, LastEdge;
static int NumbNode;

/* private functions */
static void FreeEdge (Edge edge);

/* 
 * Function: CreateNodeList
 * ------------------------
 * Allocates an empty list
 */
void CreateNodeList (int n) 
{
    LastNode = NodeList = (Node) malloc (n*sizeof (CNode));
    NumbNode = n;
}   

/* 
 * Function: FreeNodeList
 * ----------------------
 * Frees all the storage associated with the list.
 */
void FreeNodeList (void) 
{
    Node  node;
    for (node = NodeList; node != LastNode; node++) {
	free (node->x);
    }
    free (NodeList);
}

/* 
 * Function: AppendNode
 * --------------------
 * Appends a node to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.)
 */
Node AppendNode (double x, double y, double z, enum elemType type) 
{  
    Node node;

    if (LastNode == NodeList + NumbNode) return NULL; 
    node = LastNode;

    node->x = (double *) calloc (DIMENSION, sizeof(double));
    node->x[0] = x;
    node->x[1] = y;
    node->x[2] = z;
    node->type = type;

    LastNode++;
    return node;
} 


/* 
 * Function: NodeNumber
 * --------------------
 * Translates a pointer to a node into the node number in the list.
 * If node == NULL, returns the total number of nodes.
 */
int NodeNumber (Node node) 
{
    int n;

    if (node == NULL) node = LastNode;
    n = node - NodeList;
    return n;
}

/* 
 * Function: GetNode
 * -----------------
 * Returns n-th node.
 */
Node GetNode (int n)
{
    return (NodeList+n);
}

/* 
 * Function: ReadEdgeList
 * ----------------------
 * Reads an edge list from the file edgefile. 
 */
void ReadEdgeList (sf_file edgefile)
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

/* 
 * Function: WriteList
 * -----------------------
 * Writes a  node list to the file nodefile.
 * Writes an edge list to the file edgefile.
 */
void WriteList (sf_file nodefile, sf_file edgefile)
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

/* 
 * Function: NodeExec
 * ------------------
 * Recursive operations on nodes in the list.
 */
void NodeExec (void (* exec) (Node))
{
    Node node;

    for (node = NodeList; node != LastNode; node++) {
	if (node->type != EMPTY) {
	    exec (node);
	}
    } 
}

/* 
 * Function: NodeOut
 * -----------------
 * Copy a Nodes to 3 numbers
 */
void NodeOut (Node q, float *x, float *y, float *z)
{
    *x = q->x[0];
    *y = q->x[1];
    *z = q->x[2];
}

/* 
 * Function: EdgeOut
 * -----------------
 * Copy Edges to a float array. Returns the number of edges.
 */
int EdgeOut (float *e)
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

/* !!! debugging !!! */
void WriteNode (Node list)
{
    printf ("%lf %lf %lf %i\n", 
	    list->x[0], list->x[1], list->x[2], list->type);
}

void WriteEdge (Edge list)
{
    printf ("%u %u %i\n", 
	    NodeNumber (list->ends[0]), 
	    NodeNumber (list->ends[1]), list->type);
    WriteNode (list->ends[0]);
    WriteNode (list->ends[1]);
}

/* 
 * Function: CreateEdgeList
 * ------------------------
 * Allocates an empty list
 */
void CreateEdgeList (void)
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

/* 
 * Function: ClearEdgeList
 * -----------------------
 * Frees added edges in the list.
 */
void ClearEdgeList (void)
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

/* 
 * Function: FreeEdgeList
 * ----------------------
 * Frees all the storage associated with the list.
 */
void FreeEdgeList (void)
{
    FreeEdge (EdgeList);
}

/* 
 * Function: FreeEdge
 * ------------------
 * Frees all the storage associated with 
 * the edge and its descendants recursively.
 */
void FreeEdge (Edge edge)
{
    if (edge != NULL) {
	FreeEdge (edge->next);
	free (edge);
    }
}


/* 
 * Function: AppendEdge
 * --------------------
 * Appends an edge to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.)
 */ 
Edge AppendEdge (Node left, Node right, enum elemType type) 
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

/* 
 * Function: EdgeNumber
 * --------------------
 * Returns the number of edges in the list.
 */  
int EdgeNumber (void)
{
    int n;
    Edge list;

    for (n = 0, list = EdgeList; list->next != NULL; list =list->next) { 
	if (list->type != EMPTY) n++;
    }
    return n; 
}

/* 
 * Function: EdgeEnumerate
 * -----------------------
 * Give numbers to edges
 */  
void EdgeEnumerate (void)
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

/* 
 * Function: EdgeExec
 * ------------------
 * Recursive operations on edges in the list.
 */
void EdgeExec (void (* exec) (Edge))
{
    Edge edge;

    for (edge = EdgeList; edge->next != NULL; edge = edge->next) {
	if (edge->type == BOUNDARY) {
	    exec (edge);
	}
    } 
}

/*
 * Function: EdgeLength
 * --------------------
 * Returns the length of an edge
 */
double EdgeLength (Edge ab)
{
    double dx, dy;

    dx = ab->ends[1]->x[0] - ab->ends[0]->x[0];
    dy = ab->ends[1]->x[1] - ab->ends[0]->x[1];
    return (dx*dx + dy*dy);
}

/*
 * Function: MoveNode
 * --------------------
 * New position for the node
 */
void MoveNode (Node q, double x, double y)
{
    q->x[0] = x;
    q->x[1] = y;
}
