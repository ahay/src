#ifndef _celltrace_h
#define _celltrace_h

typedef struct CellTrace *celltrace;

celltrace celltrace_init (int order, int nt,
			  int nz, int nx, 
			  float dz, float dx, 
			  float z0, float x0, float* slow);
void celltrace_close (celltrace ct);
float cell_trace (celltrace ct, float* xp, float* p, int* it, float** traj);

#endif
