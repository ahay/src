/* Node operations for the tree structure. */
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

#include "node.h"

#ifndef _node_h

typedef struct CNodeCell {
    struct CNode *node;
    struct CNodeCell *link;
} *NodeCell;
/*^*/

typedef struct CNodeQueue {
    struct CNodeCell* head;
    struct CNodeCell* tail;
} *NodeQueue;
/*^*/

typedef struct CNode {
    int nparents, **parents; /* number of alive parents, immediate parents */
    int n1, n2;
    float w1, w2, t;
    struct CNodeQueue* children;
} *Node;
/*^*/

#endif

NodeQueue CreateNodeQueue (void) 
/*< start a queue >*/
{
    NodeQueue queue;
  
    queue = (NodeQueue) sf_alloc(1, sizeof(*queue));
    queue->head = NULL;
    queue->tail = NULL;

    return queue;
}

void FreeNodeQueue (NodeQueue queue) 
/*< Free allocated storage >*/
{
    NodeCell cell, next;

    for (cell = queue->head; NULL != cell; cell = next) {
	next = cell->link;
	free (cell);
    }

    free (queue);
}

Node CreateNodes (int n /* number of nodes */) 
/*< Node array allocation >*/
{
    Node nd;
    int i, j, k;

    nd = (Node) sf_alloc(n, sizeof(*nd));
    for (i=0; i < n; i++) {
	nd[i].nparents=0;
	nd[i].children = CreateNodeQueue();
	nd[i].parents = sf_intalloc2(2,2);
	for (j=0; j < 2; j++) {
	    for (k=0; k < 2; k++) {
		nd[i].parents[j][k] = -1;
	    }
	}
    }
    
    return nd;
}

void FreeNodes (Node nd, int n /* number of nodes */) 
/*< free allocated storage >*/
{
    int i;

    for (i=0; i < n; i++) {
	FreeNodeQueue(nd[i].children);
	free(nd[i].parents[0]);
	free(nd[i].parents);
    }
    free (nd);
}

void AddNode (NodeQueue queue, Node nd) 
/*< Add a node to the queue >*/
{
    NodeCell cell;

    cell = (NodeCell) sf_alloc(1,sizeof(*cell));
    cell->node = nd;
    cell->link = NULL;

    if (NULL == queue->head) {
	queue->head = cell;
    } else {
	queue->tail->link = cell;
    }
    queue->tail = cell;    
}

Node ExtractNode (NodeQueue queue) 
/*< Extract a node from the queue >*/
{
    Node nd;
    NodeCell cell;

    cell = queue->head;
    if (NULL == cell) return NULL;

    nd = cell->node;
    queue->head = cell->link;
    free (cell);

    return nd;
}


void AddChild (Node parent, int i, int j, int k, Node child) 
/*< Add a child. parent is i in the queue, [j][k] for child >*/
{
    AddNode (parent[i].children, child);
    child->parents[j][k] = i;
    child->nparents++;
}

void TraverseQueue (NodeQueue queue, void (*apply)(Node nd)) 
/*< Apply a function to every node in the queue >*/
{
    NodeCell cell;

    for (cell = queue->head; NULL != cell; cell = cell->link) {
	apply (cell->node);
    }
}

void TraverseDeleteQueue (NodeQueue queue, void (*apply)(Node nd)) 
/*< Apply a function to every node and delete the queue >*/
{
    NodeCell cell, next;

    for (cell = queue->head; NULL != cell; cell = next) {
	apply (cell->node);
	next = cell->link;
	free (cell);
    }
    queue->head = NULL;
}

/* 	$Id: node.c 1571 2005-11-21 05:50:10Z fomels $	 */
