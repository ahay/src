/*
  Copyright (C) 2012 KAUST

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

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/* MPI and OMP parameters */
extern int mpi_size;
extern int mpi_rank;
extern int omp_nth;

void exit_error(const char *format, ...)
/*< exit with error message >*/
{
	va_list args;

	fprintf(stderr,"process %4d : ",mpi_rank);

	va_start(args,format);
	vfprintf(stderr,format,args);
	va_end(args);

	fprintf(stderr,"\n");

	exit(EXIT_FAILURE);
}

void read_float(float *data,
		int size,
		const char *filename)
/*< read single precision data from file >*/
{
	FILE *file;

	if (NULL == (file = fopen(filename,"r")))
		exit_error("failed to open %s",filename);
	if (size != fread(data,sizeof(float),size,file))
		exit_error("failed reading %s",filename);
	fclose(file);
}

/* read integer data from file */
void read_int(int *data,
              int size,
              const char *filename)
{
	FILE *file;

	if (NULL == (file = fopen(filename,"r")))
		exit_error("failed to open %s",filename);
	if (size != fread(data,sizeof(int),size,file))
		exit_error("failed reading %s",filename);
	fclose(file);
}

/* fast fft size */
int fft_size(int n)
{
    int m;

    while (1) {
        m = n;
        while (!(n % 2)) n /= 2;
        while (!(n % 3)) n /= 3;
        while (!(n % 5)) n /= 5;
        while (!(n % 7)) n /= 7;
        if (1 == n) break;
        else n = m + 1;
    }
    return m;
}

void fft_expand(int n,
		int *l, /* [1] */
		int *h  /* [1] */)
/*< update boundary to fit fast fft size >*/
{
	int N,d;

	N = n + *l + *h;
	d = fft_size(N) - N;

	*l += d / 2;
	*h += d - d / 2;
}


void expand(float *q, /* [N2][N1]    */
	    const float *p,    /* [n2][n1]    */
	    const int *n,      /* {n1,n2} */
	    const int *l,      /* {l1,l2} */
	    const int *h       /* {h1,h2} */)
/*< expand 2D model >*/
{
	int i1,i2,n1,n2,l1,l2,h1,h2,N1,N2,j;
	float aa,bb;

	N1 = (n1 = n[0]) + (l1 = l[0]) + (h1 = h[0]);
	N2 = (n2 = n[1]) + (l2 = l[1]) + (h2 = h[1]);

	for (i2=0; i2 < n2; i2++) { /* extend axis 1 */
		aa = p[i2*n1];
		bb = p[i2*n1 + n1-1];
		for (i1=0; i1 < N1; i1++) {
			j = (l2+i2) * N1 + i1;
			if      (i1 < l1) 
				q[j] = aa;
			else if (i1 < l1 + n1)
				q[j] = p[i2*n1 + i1-l1];
			else
				q[j] = bb;
		}
	}

	for (i1=0; i1 < N1; i1++) { /* extend axis 2 */
		for (aa = q[l2*N1 + i1], i2=0; i2 < l2; i2++)
			q[i2*N1 + i1] = aa;
		for (bb = q[(l2+n2-1)*N1 + i1], i2=l2+n2; i2 < N2; i2++)
			q[i2*N1 + i1] = bb;
	}
}


void compute_k(float *k, /* [nk] */
	       int nk)
/*< compute wavenumbers
   0     to N/2 positive
   N/2+1 to N-1 negative >*/
{
	int i;
	float dk;

	for (k[0]=0, dk=1./nk, i=1; i < nk; i++) {
		k[i] = k[i-1] + dk;
		if (nk/2+1 == i) k[i] -= 1.; 
	}
}

int cut_offset(float *data,   /* [nr][n4] NULL if cut offset only */
               float *cr,     /* [nr][2 ] */
               const int *n,  /*  {nr,n4} */
               const float *x /*  {x0,xe} */)
/*< cut out-of-domain offset and traces >*/
{
    int i,j;
    int nr = n[0];
    int n4 = n[1];
    float x0 = x[0];
    float xe = x[1];

    if (cr[1] < cr[nr*2-1]) {
        i = 0;    while (cr[i*2+1] < x0) i++;
        j = nr-1; while (cr[j*2+1] > xe) j--;
        printf("process %3d : i=%d j=%d skip first %d and last %d traces\n",mpi_rank,i,j,i,nr-1-j);
        cr = &cr[i*2];
        if (NULL != data) data = &data[i*n4];
        nr = j - i + 1;
    } else {
        i = nr-1; while (cr[i*2+1] < x0) i--;
        j = 0;    while (cr[j*2+1] > xe) j++;
        printf("process %3d : i=%d j=%d SKIP first %d and last %d traces\n",mpi_rank,i,j,j,nr-1-i);
        cr = &cr[j*2];
        if (NULL != data) data = &data[j*n4];
        nr = i - j + 1;
    }
    
    return nr;
}
