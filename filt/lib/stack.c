#include <stdio.h>

#include "c99.h"

#include "stack.h"
#include "alloc.h"
#include "error.h"

struct entry
{
    void* data;
    int type;
};

struct sf_Stack {
    int size;
    struct entry *entry, *top;
};

sf_stack sf_stack_init (int size)
{
    sf_stack s;

    s = (sf_stack) sf_alloc (1, sizeof(*s));
    s->size = size;
    s->entry = (struct entry*) sf_alloc(size,sizeof(struct entry));
    s->top = s->entry-1;

    return s;
}

void sf_stack_print (sf_stack s)
{
    struct entry *e;

    fprintf(stderr,"----------------------\n");
    for (e=s->entry; e <= s->top; e++) {
	fprintf(stderr,"type=%d data=%p\n",e->type,e->data);
    }
    fprintf(stderr,"----------------------\n");   
}

int sf_stack_get (sf_stack s) {
    return (s->top - s->entry);
}

void sf_stack_set (sf_stack s, int pos) {
    s->top = s->entry + pos;
}

/* requires unique data for each push */
void sf_push(sf_stack s, void *data, int type)
{
    if (s->top == s->entry + s->size) sf_error("%s: stack is full",__FILE__); 

    s->top++;
    s->top->type = type;
    s->top->data = data;
}

void* sf_pop(sf_stack s)
{
    struct entry *old;

    old = s->top;
    s->top--;

    return old->data;
}

bool sf_full (sf_stack s)
{
    return (s->top >= s->entry);
}

int sf_top(sf_stack s)
{
    return s->top->type;
}

void sf_stack_close(sf_stack s)
{
    free (s->entry);
    free (s);
}

/* 	$Id: stack.c,v 1.2 2003/09/29 14:34:56 fomels Exp $	 */
