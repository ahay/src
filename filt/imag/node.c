#include <stdlib.h>

#include <rsf.h>

#include "node.h"

NodeQueue CreateNodeQueue (void) {
    NodeQueue queue;
  
    queue = (NodeQueue) sf_alloc(1, sizeof(*queue));
    queue->head = NULL;
    queue->tail = NULL;

    return queue;
}

void FreeNodeQueue (NodeQueue queue) {
    NodeCell cell, next;

    for (cell = queue->head; NULL != cell; cell = next) {
	next = cell->link;
	free (cell);
    }

    free (queue);
}

Node CreateNodes (int n, int order) {
    Node nd;
    int i, j, k;

    nd = (Node) sf_alloc(n, sizeof(*nd));
    for (i=0; i < n; i++) {
	nd[i].nparents=0;
	nd[i].children = CreateNodeQueue();
	nd[i].parents = sf_intalloc2(order,order);
	for (j=0; j < order; j++) {
	    for (k=0; k < order; k++) {
		nd[i].parents[j][k] = -1;
	    }
	}
    }
    
    return nd;
}

void FreeNodes (Node nd, int n) {
    int i;

    for (i=0; i < n; i++) {
	FreeNodeQueue(nd[i].children);
	free(nd[i].parents[0]);
	free(nd[i].parents);
    }
    free (nd);
}

void AddNode (NodeQueue queue, Node nd) {
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

Node ExtractNode (NodeQueue queue) {
    Node nd;
    NodeCell cell;

    cell = queue->head;
    if (NULL == cell) return NULL;

    nd = cell->node;
    queue->head = cell->link;
    free (cell);

    return nd;
}

/* parent is i in the queue, [j][k] for child */
void AddChild (Node parent, int i, int j, int k, Node child) {
    AddNode (parent[i].children, child);
    child->parents[j][k] = i;
    child->nparents++;
}

void TraverseQueue (NodeQueue queue, void (*apply)(Node nd)) {
    NodeCell cell;

    for (cell = queue->head; NULL != cell; cell = cell->link) {
	apply (cell->node);
    }
}

void TraverseDeleteQueue (NodeQueue queue, void (*apply)(Node nd)) {
    NodeCell cell, next;

    for (cell = queue->head; NULL != cell; cell = next) {
	apply (cell->node);
	next = cell->link;
	free (cell);
    }
    queue->head = NULL;
}
