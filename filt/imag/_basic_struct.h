/*
 * File: basic_struct.h
 * --------------------
 * Contains structure declaration for "Node", "Edge", and "Triangle". 

 * Course CS368, 05/19/97, 
 * Project --- Conforming Delauney Triangulation.
 */
#ifndef _basic_struct_h
#define _basic_struct_h

typedef struct CNode {
    double *x;
    enum elemType type;
} CNode;


typedef struct CEdge {
    struct CNode *ends[2];
    struct CTriangle* face[2];
    enum elemType type;
    int num;
    struct CEdge *next;
} CEdge;

typedef struct CTriangle {
    struct CEdge *edge[3];
    struct CTriangle **child;
    enum elemType clone;
} *Triangle;

#endif




