/*
 * Interface for node and edge list manipulations
 */
#ifndef _list_struct_h
#define _list_struct_h

/* number of dimensions */
#define DIMENSION 3

/* 
 * Type: elemType
 * --------------
 * Flag to distinguish elements
 */
enum elemType {EMPTY = -1,BOUNDARY,ADDED};

/* 
 * Type: Node
 * ----------
 * Abstract data type for Nodes
 */
typedef struct CNode *Node;

/* 
 * Type: Edge
 * ----------
 * Abstract data type for Edges
 */
typedef struct CEdge *Edge;

/* 
 * Function: CreateNodeList
 * ------------------------
 * Allocates an empty list
 */
void CreateNodeList (int n);

/* 
 * Function: FreeNodeList
 * ----------------------
 * Frees all the storage associated with the list.
 */
void FreeNodeList (void);

/* 
 * Function: AppendNode
 * --------------------
 * Appends a node to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.)
 */
Node AppendNode (double x, double y, double z, enum elemType type);

/* 
 * Function: ReadEdgeList
 * ----------------------
 * Reads an edge list from the file edgefile. 
 * If filename == "in", reads from stdin. 
 * edgefile should have the following form:
 * 
 * number_of_edges
 * i_0 j_0 
 * i_1 j_1 
 * ...
 *
 * i's and j's correspond to the node number in the list.
 */
void ReadEdgeList  (sf_file edgefile);

/* 
 * Function: WriteList
 * -----------------------
 * Writes a  node list to the file nodefile.
 * Writes an edge list to the file edgefile.
 * If filename == "out", writes to stdout. 
 * File formats defined in ReadList.
 */
void WriteList (sf_file nodefile, sf_file edgefile);

/* 
 * Function: CreateEdgeList
 * ------------------------
 * Allocates an empty list
 */
void CreateEdgeList (void);

/* 
 * Function: FreeEdgeList
 * ----------------------
 * Frees all the storage associated with the list.
 */
void FreeEdgeList (void);

/* 
 * Function: ClearEdgeList
 * -----------------------
 * Frees added edges in the list.
 */
void ClearEdgeList (void);

/* 
 * Function: AppendEdge
 * --------------------
 * Appends an edge to the end of the list. 
 * Returns a pointer to the inserted node 
 * (useful for incremental insertions.)
 */
Edge AppendEdge (Node left, Node right, enum elemType type);

/* 
 * Function: NodeNumber
 * --------------------
 * Translates a pointer to a node into the node number in the list.
 * If node == NULL, returns the total number of nodes.
 */
int NodeNumber (Node node);

/* 
 * Function: NodeExec
 * ------------------
 * Recursive operations on nodes in the list.
 */
void NodeExec (void (* exec) (Node));

/* 
 * Function: NodeExec
 * ------------------
 * Recursive operations on nodes in the list.
 */
void EdgeExec (void (* exec) (Edge));

/* 
 * Function: EdgeNumber
 * --------------------
 * Returns the number of edges in the list.
 */  
int EdgeNumber (void);


/* 
 * Function: NodeOut
 * -----------------
 * Copy a Node to 3 numbers.
 */
void NodeOut (Node q, float *x, float *y, float *z);

/* 
 * Function: EdgeOut
 * -----------------
 * Copy Edges to a float array. Returns the number of edges.
 */
int EdgeOut (float *e);

/*
 * Function: EdgeLength
 * --------------------
 * Returns the length of an edge
 */
double EdgeLength (Edge ab);

/* 
 * Function: GetNode
 * -----------------
 * Returns n-th node.
 */
Node GetNode (int n);

/*
 * Function: MoveNode
 * --------------------
 * New position for the node
 */
void MoveNode (Node q, double x, double y); 

/* 
 * Function: EdgeEnumerate
 * -----------------------
 * Give numbers to edges
 */  
void EdgeEnumerate (void);

#endif
