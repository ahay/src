#ifndef _pfafft_h
#define _pfafft_h

#include "c99.h"

#ifndef __cplusplus

/* Prime Factor FFTs */
int sf_npfa (int nmin);
int sf_npfao (int nmin, int nmax);
int sf_npfar (int nmin);
int sf_npfaro (int nmin, int nmax);
void sf_pfacc (int isign, int n, float complex z[]);
void sf_pfarc (int isign, int n, float rz[], float complex cz[]);
void sf_pfacr (int isign, int n, float complex cz[], float rz[]);
void sf_pfa2cc (int isign, int idim, int n1, int n2, float complex z[]);
void sf_pfa2rc (int isign, int idim, int n1, int n2, 
		float rz[], float complex cz[]);
void sf_pfa2cr (int isign, int idim, int n1, int n2, 
		float complex cz[], float rz[]);
void sf_pfamcc (int isign, int n, int nt, int k, int kt, float complex z[]);

#endif 

#endif

/* 	$Id$	 */

