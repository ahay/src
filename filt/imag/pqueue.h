#ifndef _pqueue_h
#define _pqueue_h

enum {FMM_IN, FMM_FRONT, FMM_OUT};

/* Initialize heap with the maximum size */
void pqueue_init (int n);

/* Set starting values */
void pqueue_start (void);

/* Free the allocated storage */
void pqueue_close (void);

/* Insert an element */
void pqueue_insert (float* v);

/* Extract the time element */
float* pqueue_extract (void);

#endif
