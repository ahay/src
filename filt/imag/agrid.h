#ifndef _agrid_h
#define _agrid_h

typedef struct AGrid *agrid;

agrid agrid_init (int n, int nd, int max);
void agrid_set (agrid grid, float** dat);
void agrid_close (agrid grid);
void agrid_interp (agrid grid, int i, float x, float* f);
void fill_grid (agrid grid, float min1, float max1, float min2, float max2,
		void* dat, void (*fill)(float x, void* dat, float* f));
int grid_size (agrid grid);
float** write_grid (agrid grid);


#endif
