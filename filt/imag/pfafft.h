#ifndef _pfafft_h
#define _pfafft_h

#include <rsf.h>

#ifndef PI
#define PI (3.141592653589793)
#endif

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, float complex z[]);
void pfarc (int isign, int n, float rz[], float complex cz[]);
void pfacr (int isign, int n, float complex cz[], float rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, float complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, 
	     float rz[], float complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, 
	     float complex cz[], float rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, float complex z[]);

#endif 
