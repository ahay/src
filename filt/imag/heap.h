/* 
 * Copyright (c) 1997 Stanford Exploration Project
 * All Rights Reserved
 *
 * File: heap.h 
 */

/*
 * This file defines the interface for the
 * simple heap priority queue
 *
 * Author: Sergey Fomel
 */
#ifndef _heap_h
#define _heap_h

typedef struct Point {
  double v;
  void* d;
  struct Point **h;
} Point; 

/* initialize it with the maximum number of elements */
void heap_init (int nx);

/* free the storage */
void heap_close (void);

/* insert a point */
void heap_insert (Point* v);

/* extract the maximum */
Point* heap_extract (void);

/* restore the heap the value has been altered */
void heap_update (Point* v);

#endif
