#ifndef _node_h
#define _node_h

typedef struct CNodeCell {
    struct CNode *node;
    struct CNodeCell *link;
} *NodeCell;

typedef struct CNodeQueue {
    struct CNodeCell* head;
    struct CNodeCell* tail;
} *NodeQueue;
 
typedef struct CNode {
    int nparents, **parents; /* number of alive parents, immediate parents */
    int n1, n2;
    float w1, w2, t;
    struct CNodeQueue* children;
} *Node;

NodeQueue CreateNodeQueue (void);

void FreeNodeQueue (NodeQueue queue);

Node CreateNodes (int n, int order);

void FreeNodes (Node nd, int n);

void AddNode (NodeQueue queue, Node nd);

void AddChild (Node parent, int i, int j, int k, Node child);

void TraverseQueue (NodeQueue queue, void (*apply)(Node nd));

void TraverseDeleteQueue (NodeQueue queue, void (*apply)(Node nd));

Node ExtractNode (NodeQueue queue);

#endif
