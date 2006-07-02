/* Generic stack (FILO) structure operations. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <stdio.h>

#include "_bool.h"
/*^*/

#include "stack.h"
#include "alloc.h"
#include "error.h"

#ifndef _sf_stack_h

typedef struct sf_Stack *sf_stack;
/*^*/

#endif

struct entry
{
    void* data;
    int type;
};

struct sf_Stack {
    size_t size;
    struct entry *entry, *top;
};

sf_stack sf_stack_init (size_t size)
/*< create a stack >*/
{
    sf_stack s;

    s = (sf_stack) sf_alloc (1, sizeof(*s));
    s->size = size;
    s->entry = (struct entry*) sf_alloc(size,sizeof(struct entry));
    s->top = s->entry-1;

    return s;
}

void sf_stack_print (sf_stack s)
/*< print out a stack (for debugging) >*/ 
{
    struct entry *e;

    fprintf(stderr,"----------------------\n");
    for (e=s->entry; e <= s->top; e++) {
	fprintf(stderr,"type=%d data=%p\n",e->type,e->data);
    }
    fprintf(stderr,"----------------------\n");   
}

int sf_stack_get (sf_stack s) 
/*< extract stack length >*/
{
    return (s->top - s->entry);
}

void sf_stack_set (sf_stack s, int pos) 
/*< set stack position >*/
{
    s->top = s->entry + pos;
}

void sf_push(sf_stack s, void *data, int type)
/*< push data into stack (requires unique data for each push) >*/
{
    if (s->top == s->entry + s->size) sf_error("%s: stack is full",__FILE__); 

    s->top++;
    s->top->type = type;
    s->top->data = data;
}

void* sf_pop(sf_stack s)
/*< pop data from stack >*/ 
{
    struct entry *old;

    old = s->top;
    s->top--;

    return old->data;
}

bool sf_full (sf_stack s)
/*< test if the stack is full >*/
{
    return (bool) (s->top >= s->entry);
}

int sf_top(sf_stack s)
/*< return the top type >*/
{
    return s->top->type;
}

void sf_stack_close(sf_stack s)
/*< free allocated memory >*/
{
    free (s->entry);
    free (s);
}

/* 	$Id$	 */
