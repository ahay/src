#ifndef _sf_stack_h
#define _sf_stack_h

#include "c99.h"

typedef struct sf_Stack *sf_stack;

sf_stack sf_stack_init (int size);
int sf_stack_get (sf_stack s);
void sf_stack_set (sf_stack s, int pos);
void sf_push(sf_stack s, void *data, int type);
void* sf_pop(sf_stack s);
int sf_top(sf_stack s);
bool sf_full (sf_stack s);
void sf_stack_close(sf_stack s);
void sf_stack_print (sf_stack s);

#endif

/* 	$Id: stack.h,v 1.2 2003/09/29 14:34:56 fomels Exp $	 */
