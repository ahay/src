#ifndef _tree_h
#define _tree_h

void tree_init (int order,
		int nz1, int nx1, int na1, int nt1, 
		float dz1, float dx1, float dp1, 
		float z01, float x01, float p01, 
		float** vel, float** val);
void tree_traverse(void);
void tree_build(int method);
void tree_close(void);
void tree_print (void);

#endif
