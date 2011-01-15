/* Linked list for use in conjugate-direction-type methods. */
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

#include "llist.h"
#include "alloc.h"

#include <stdio.h>


#ifndef _sf_llist_h

typedef struct sf_List *sf_list;
/* abstract data type */
/*^*/

#endif

struct Entry {
    float *g;
    double gn;
    struct Entry *next;
};

struct sf_List {
    struct Entry *first, *current;
    int depth;
};

sf_list sf_llist_init(void)
/*< create an empty list >*/
{
    sf_list l;

    l = (sf_list) sf_alloc(1,sizeof(*l));
    l->first = (struct Entry *) sf_alloc(1,sizeof(struct Entry));
    l->depth = 0;
    sf_llist_rewind(l);

    return l;
}

void sf_llist_rewind(sf_list l)
/*< return to the start >*/
{
    l->current = l->first;
}

int sf_llist_depth(sf_list l)
/*< return list depth >*/
{
    return l->depth;
}

void sf_llist_add(sf_list l, float *g, double gn) 
/*< add an entry in the list >*/
{    
    l->current->g = g;
    l->current->gn = gn;
    l->current->next = (struct Entry *) sf_alloc(1,sizeof(struct Entry)); /* never deallocated */
    l->depth++;
}

void sf_llist_down(sf_list l, float **g, double *gn)
/*< extract and entry from the list >*/
{
    *g = l->current->g;
    *gn = l->current->gn;
    l->current = l->current->next;
}

void sf_llist_close(sf_list l) 
/*< free allocated storage >*/
{
    int i,depth;

    depth = l->depth;

    for (i=0; i < depth; i++) {
	sf_llist_chop(l);
    }

    free(l->first);
    free(l);
}

void sf_llist_chop(sf_list l) 
/*< free the top entry from the list >*/
{
    sf_llist_rewind(l); 

    l->first = l->current->next;
    if (NULL != l->current->g)
    {
    	free (l->current->g);
    	l->current->g = NULL;
    }
    free (l->current);
    l->depth--;

}

    
