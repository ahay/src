#ifndef _acell_h
#define _acell_h

typedef struct ACell *acell;

acell acell_init (void);
void acell_set (acell cell, float* dat1, float* dat2);
void free_acell (acell cell);
void acells_init (int n, int max);
void free_acells (void);
void interp_acell (acell cell, float x, float* f);
void fill_cell (acell cell, int i, void* dat, 
		void (*fill)(float x, void* dat, float* f));
void split_cell (acell cell, 
		 float min1, float max1, 
		 float min2, float max2, int i, void* dat,
		 void (*fill)(float x, void* dat, float* f));
void fill_node (acell cell, int i, float x, void* dat,
		void (*fill)(float x, void* dat, float* f));
int cell_size (acell cell);
float* get_node (acell cell, int i);
void flat_cell (acell cell, acell *flat, float** dat);

#endif
