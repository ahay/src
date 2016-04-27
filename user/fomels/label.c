/* Structure for labeling algorithm */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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
#include <rsf.h>

#include "label.h"

static int *table;
static int label;
static int n;

static int find_root(int i);
static void set_root(int i, int r);

void label_init(int maxl /* maximum number of labels */)
/*< initialize >*/
{
    table = sf_intalloc(maxl);
    label = 0;
    n = 0;
}

void label_close(void)
/*< free allocated storage >*/
{
    free (table);
}

static int find_root(int i)
{
    while (table[i] < i) 
	i = table[i];
    return i;
}

static void set_root(int i, int r)
{
    int j;
    
    while (table[i] < i) {
	j = table[i];
	table[i] = r;
	i = j;
    }
    table[i] = r;
}

void label_union(int i, int j)
/*< join two labels >*/
{
    int ri, rj;
    
    if (i != j) {
	ri = find_root(i);
	rj = find_root(j);
	if (ri > rj) ri = rj;
	set_root(i,ri);
	set_root(j,ri);
    }
}

int label_new(void)
/*< create new label >*/
{
    int r;
    
    r = label;
    label++;

    table[n] = r;
    n++;
    
    return r;
}

void label_flatten(void)
/*< flatten label list >*/
{
    int i;

    for (i=1; i < n; i++) {
	table[i] = table[table[i]];
    }
}

int label_find(int i)
/*< find label >*/
{
    int r;

    r = find_root(i);
    set_root(i,r);
    return r;
}
    
  
