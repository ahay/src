#ifndef _node_h
#define _node_h

typedef struct CNodeList {
  int nitems, ntotal;
  struct CNode** list;
} *NodeList;
 
typedef struct CNode {
  int nparents, **parents; /* number of alive parents, immediate parents */
  float w1, w2, t;
  struct CNodeList* children;
} *Node;

NodeList CreateNodeList (int n);

Node CreateNodes (int n, int order);

void FreeNodes (Node nd, int n);

void AddNode (NodeList list, Node nd);

void AddChild (Node parent, int i, int j, int k, Node child);

#endif
