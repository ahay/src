#include <stdlib.h>

#include <rsf.h>

#include "node.h"

NodeList CreateNodeList (int n) {
    NodeList list;
  
    list = (NodeList) sf_alloc(1, sizeof(*list));
    list->nitems = 0;
    list->ntotal = n;
    list->list = (Node*) sf_alloc(n,sizeof(Node));

    return list;
}

Node CreateNodes (int n) {
    Node nd;
    int i;

    nd = (Node) sf_alloc(n, sizeof(*nd));
    for (i=0; i < n; i++) {
	nd[i].nparents=0;
	nd[i].children = CreateNodeList(1);
	nd[i].parents[0] = -1;
	nd[i].parents[1] = -1;
	nd[i].parents[2] = -1;
	nd[i].parents[3] = -1;
    }

    return nd;
}

void FreeNodes (Node nd, int n) {
    int i;

    for (i=0; i < n; i++) {
	free(nd[i].children->list);
	free(nd[i].children);
    }
    free (nd);
}

void AddChild (Node parent, int i, int k, Node child) {
    AddNode (parent[i].children, child);
    child->parents[k] = i;
    child->nparents++;
}

void AddNode (NodeList list, Node nd) {
    list->list[list->nitems] = nd;
    list->nitems++;
    if (list->nitems == list->ntotal) {
	list->ntotal *= 2;
	list->list = (Node*) sf_realloc (list->list,list->ntotal,sizeof(Node));
    }
}

