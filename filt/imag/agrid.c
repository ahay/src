#include <rsf.h>

#include "agrid.h"
#include "acell.h"

struct AGrid {
    int n;
    acell *cells;
    int maxchild;
};

agrid agrid_init (int n, int nd, int max)
{
    int i;
    agrid grid;

    grid = (agrid) sf_alloc (1,sizeof(*grid));

    acells_init (nd,max);

    grid->cells = (acell*) sf_alloc (n-1,sizeof(acell));
    grid->n = n;

    for (i=0; i < n-1; i++) {
	grid->cells[i] = acell_init ();
    }

    return grid;
}

void agrid_set (agrid grid, float** dat) 
{
    int i;

    for (i=0; i < grid->n-1; i++) {
	acell_set (grid->cells[i],dat[i],dat[i+1]);
    }
}

void agrid_close (agrid grid)
{
    int i;

    for (i=0; i < grid->n-1; i++) {
	free_acell (grid->cells[i]);
    }
    free (grid->cells);
/*    free_acells(); */
    free (grid);
}

void agrid_interp (agrid grid, int i, float x, float* f) {
    interp_acell (grid->cells[i],x, f);
}

/* x is relative to the start of the grid */
void fill_grid (agrid grid, float min1, float max1, float min2, float max2,
		void* dat, void (*fill)(float x, void* dat, float* f)) {
    int i, n;

    n = grid->n-1;
    for (i=0; i < n; i++) {
	fill_cell (grid->cells[i], i, dat, fill);
    }
    fill_node (grid->cells[n-1], 1, (float) n, dat, fill);
    for (i=0; i < n; i++) {
	split_cell (grid->cells[i], min1, max1, min2, max2, i, dat, fill);
    }
}

int grid_size (agrid grid) {
    int i, n, size;

    n = grid->n-1;
    size = 0;
    for (i=0; i < n; i++) {
	size += cell_size (grid->cells[i]);
    }
    size++;

    return size;
}

float** write_grid (agrid grid) {
    acell* flat;
    float** dat;
    int i, k, n, *size, maxsize, totsize;

    n = grid->n-1;

    size = sf_intalloc (n);

    totsize = 0;
    maxsize = 0;
    for (i=0; i < n; i++) {
	size[i] = cell_size (grid->cells[i]);
	if (size[i] > maxsize) maxsize = size[i];
	totsize += size[i];
    }
    totsize++;

    dat = (float**) sf_alloc(totsize,sizeof(float*));
    flat = (acell*) sf_alloc(maxsize,sizeof(acell));
    for (i=0, k=0; i < n; k+= size[i], i++) {
	flat_cell (grid->cells[i],flat,dat+k);
    }
    dat[k] = get_node(grid->cells[n-1],1);
    free (size);
    free (flat);

    return dat;
}

/* 	$Id: agrid.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */

