#ifndef _node_h
#define _node_h

typedef struct CNodeList {
  int nitems, ntotal;
  struct CNode** list;
} *NodeList;
 
typedef struct CNode {
  int nparents;
  int parents[4];
  float w1, w2, t;
  struct CNodeList* children;
} *Node;

NodeList CreateNodeList (int n);

Node CreateNodes (int n);

void FreeNodes (Node nd, int n);

void AddNode (NodeList list, Node nd);

void AddChild (Node parent, int i, int k, Node child);

#endif
