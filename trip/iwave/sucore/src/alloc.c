/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
ALLOC - Allocate and free multi-dimensional arrays

alloc1		allocate a 1-d array
realloc1	re-allocate a 1-d array
free1		free a 1-d array
alloc2		allocate a 2-d array
free2		free a 2-d array
alloc3		allocate a 3-d array
free3		free a 3-d array
alloc4		allocate a 4-d array
free4		free a 4-d array
alloc5		allocate a 5-d array
free5		free a 5-d array
alloc6		allocate a 6-d array
free6		free a 6-d arrayalloc1int	
allocate a 1-d array of ints
realloc1int	re-allocate a 1-d array of ints
free1int	free a 1-d array of ints
alloc2int	allocate a 2-d array of ints
free2int	free a 2-d array of ints
alloc3int	allocate a 3-d array of ints
free3int	free a 3-d array of ints
alloc1float	allocate a 1-d array of floats
realloc1float	re-allocate a 1-d array of floats
free1float	free a 1-d array of floats
alloc2float	allocate a 2-d array of floats
free2float	free a 2-d array of floats
alloc3float	allocate a 3-d array of floats
free3float	free a 3-d array of floats
alloc4float	allocate a 4-d array of floats 
free4float      free a 4-d array of floats 
alloc5float     allocate a 5-d array of floats 
free5float      free a 5-d array of floats 
alloc6float     allocate a 6-d array of floats 
free6float      free a 6-d array of floats 
alloc4int       allocate a 4-d array of ints 
free4int        free a 4-d array of ints 
alloc5int       allocate a 5-d array of ints 
free5int        free a 5-d array of ints 
alloc5uchar	allocate a 5-d array of unsigned chars 
free5uchar	free a 5-d array of unsiged chars 
alloc5ushort    allocate a 5-d array of unsigned shorts 
free5ushort     free a 5-d array of unsiged shorts
alloc6ushort    allocate a 6-d array of unsigned shorts 
free6ushort     free a 6-d array of unsiged shorts
alloc1double	allocate a 1-d array of doubles
realloc1double	re-allocate a 1-d array of doubles
free1double	free a 1-d array of doubles
alloc2double	allocate a 2-d array of doubles
free2double	free a 2-d array of doubles
alloc3double	allocate a 3-d array of doubles
free3double	free a 3-d array of doubles
alloc1complex	allocate a 1-d array of complexs
realloc1complex	re-allocate a 1-d array of complexs
free1complex	free a 1-d array of complexs
alloc2complex	allocate a 2-d array of complexs
free2complex	free a 2-d array of complexs
alloc3complex	allocate a 3-d array of complexs
free3complex	free a 3-d array of complexs

alloc1dcomplex   allocate a 1-d array of complexs
realloc1dcomplex re-allocate a 1-d array of complexs
free1dcomplex    free a 1-d array of complexs
alloc2dcomplex   allocate a 2-d array of complexs 
free2dcomplex    free a 2-d array of complexs
alloc3dcomplex   allocate a 3-d array of complexs  
free3dcomplex    free a 3-d array of complexs

******************************************************************************
Function Prototypes:
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void free1 (void *p);
void **alloc2 (size_t n1, size_t n2, size_t size);
void free2 (void **p);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void free3 (void ***p);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void free4 (void ****p);
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void free5 (void *****p);
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6, 
                  size_t size);
void free6 (void ******p);
int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
void free1int (int *p);
int **alloc2int (size_t n1, size_t n2);
void free2int (int **p);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
void free3int (int ***p);
float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
void free1float (float *p);
float **alloc2float (size_t n1, size_t n2);
void free2float (float **p);
float ***alloc3float (size_t n1, size_t n2, size_t n3);
void free3float (float ***p);
float ****alloc4float (size_t n1, size_t n2, size_t n3, size_t n4);
void free4float (float ****p);
float *****alloc5float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5float (float *****p);
float ******alloc6float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, 
                         size_t n6);
void free6float (float ******p);
int ****alloc4int (size_t n1, size_t n2, size_t n3, size_t n4);
void free4int (int ****p);
int *****alloc5int (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5int (int *****p);
unsigned char *****alloc5uchar(size_t n1,size_t n2,size_t n3,size_t n4,
	size_t n5);
void free5uchar(unsigned char *****p);
unsigned short *****alloc5ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
void free5ushort(unsigned short *****p);
unsigned short ******alloc6ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5,size_t n6);
void free6ushort(unsigned short ******p);
double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
void free1double (double *p);
double **alloc2double (size_t n1, size_t n2);
void free2double (double **p);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
void free3double (double ***p);
complex *alloc1complex (size_t n1);
complex *realloc1complex (complex *v, size_t n1);
void free1complex (complex *p);
complex **alloc2complex (size_t n1, size_t n2);
void free2complex (complex **p);
complex ***alloc3complex (size_t n1, size_t n2, size_t n3);
void free3complex (complex ***p);

complex *alloc1dcomplex (size_t n1);
complex *realloc1dcomplex (dcomplex *v, size_t n1);
void free1dcomplex (dcomplex *p);
complex **alloc2dcomplex (size_t n1, size_t n2);
void free2dcomplex (dcomplex **p);
complex ***alloc3dcomplex (size_t n1, size_t n2, size_t n3);
void free3dcomplex (dcomplex ***p);


******************************************************************************
Notes:
The functions defined below are intended to simplify manipulation
of multi-dimensional arrays in scientific programming in C.  These
functions are useful only because true multi-dimensional arrays
in C cannot have variable dimensions (as in FORTRAN).  For example,
the following function IS NOT valid in C:
	void badFunc(a,n1,n2)
	float a[n2][n1];
	{
		a[n2-1][n1-1] = 1.0;
	}
However, the following function IS valid in C:
	void goodFunc(a,n1,n2)
	float **a;
	{
		a[n2-1][n1-1] = 1.0;
	}
Therefore, the functions defined below do not allocate true
multi-dimensional arrays, as described in the C specification.
Instead, they allocate and initialize pointers (and pointers to 
pointers) so that, for example, a[i2][i1] behaves like a 2-D array.

The array dimensions are numbered, which makes it easy to add 
functions for arrays of higher dimensions.  In particular,
the 1st dimension of length n1 is always the fastest dimension,
the 2nd dimension of length n2 is the next fastest dimension,
and so on.  Note that the 1st (fastest) dimension n1 is the 
first argument to the allocation functions defined below, but 
that the 1st dimension is the last subscript in a[i2][i1].
(This is another important difference between C and Fortran.)

The allocation of pointers to pointers implies that more storage
is required than is necessary to hold a true multi-dimensional array.
The fraction of the total storage allocated that is used to hold 
pointers is approximately 1/(n1+1).  This extra storage is unlikely
to represent a significant waste for large n1.

The functions defined below are significantly different from similar 
functions described by Press et al, 1988, NR in C.
In particular, the functions defined below:
	(1) Allocate arrays of arbitrary size elements.
	(2) Allocate contiguous storage for arrays.
	(3) Return NULL if allocation fails (just like malloc).
	(4) Do not provide arbitrary lower and upper bounds for arrays.

Contiguous storage enables an allocated multi-dimensional array to
be passed to a C function that expects a one-dimensional array.
For example, to allocate and zero an n1 by n2 two-dimensional array
of floats, one could use
	a = alloc2(n1,n2,sizeof(float));
	zeroFloatArray(n1*n2,a[0]);
where zeroFloatArray is a function defined as
	void zeroFloatArray(int n, float *a)
	{
		int i;
		for (i=0; i<n; i++)
			a[i] = 0.0;
	}

Internal error handling and arbitrary array bounds, if desired,
should be implemented in functions that call the functions defined 
below, with the understanding that these enhancements may limit 
portability.

******************************************************************************
Author:    	Dave Hale, Colorado School of Mines, 12/31/89
                Zhaobo Meng, added 4D, 5D and 6D functions, 1996
*****************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"

/* allocate a 1-d array */
void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	size_t i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 4-d array */
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
	size_t i4,i3,i2;
	void ****p;

	if ((p=(void****)malloc(n4*sizeof(void***)))==NULL)
		return NULL;
	if ((p[0]=(void***)malloc(n4*n3*sizeof(void**)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void**)malloc(n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void*)malloc(n4*n3*n2*n1*size))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i4=0; i4<n4; i4++) {
		p[i4] = p[0]+i4*n3;
		for (i3=0; i3<n3; i3++) {
			p[i4][i3] = p[0][0]+n2*(i3+n3*i4);
			for (i2=0; i2<n2; i2++)
				p[i4][i3][i2] = (char*)p[0][0][0]+
						size*n1*(i2+n2*(i3+n3*i4));
		}
	}
	return p;
}

/* free a 4-d array */
void free4 (void ****p)
{
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 5-d array */
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
	size_t i5,i4,i3,i2;
	void *****p;

	if ((p=(void*****)malloc(n5*sizeof(void****)))==NULL)
		return NULL;
	if ((p[0]=(void****)malloc(n5*n4*sizeof(void***)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void***)malloc(n5*n4*n3*sizeof(void**)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void**)malloc(n5*n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0]=(void*)malloc(n5*n4*n3*n2*n1*size))==NULL) {
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i5=0; i5<n5; i5++) {
		p[i5] = p[0]+n4*i5;
		for (i4=0; i4<n4; i4++) {
			p[i5][i4] = p[0][0]+n3*(i4+n4*i5);
			for (i3=0; i3<n3; i3++) {
				p[i5][i4][i3] = p[0][0][0]+n2*(i3+n3*(i4+n4*i5));
				for (i2=0; i2<n2; i2++)
					p[i5][i4][i3][i2] = (char*)p[0][0][0][0]+
						size*n1*(i2+n2*(i3+n3*(i4+n4*i5)));
			}
		}
	}
	return p;
}

/* free a 5-d array */
void free5 (void *****p)
{
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 6-d array */
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size)
{
	size_t i6,i5,i4,i3,i2;
	void ******p;

	if ((p=(void******)malloc(n6*sizeof(void*****)))==NULL)
		return NULL;

	if ((p[0]=(void*****)malloc(n6*n5*sizeof(void****)))==NULL) {
                free(p);
		return NULL;
        }

	if ((p[0][0]=(void****)malloc(n6*n5*n4*sizeof(void***)))==NULL) {
		free(p[0]);
                free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void***)malloc(n6*n5*n4*n3*sizeof(void**)))==NULL) {
		free(p[0][0]);
                free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0]=(void**)malloc(n6*n5*n4*n3*n2*sizeof(void*)))==NULL) {
	        free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0][0]=(void*)malloc(n6*n5*n4*n3*n2*n1*size))==NULL) {
	        free(p[0][0][0][0]);
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}

        for (i6=0; i6<n6; i6++) {
                p[i6] = p[0]+n5*i6;
	        for (i5=0; i5<n5; i5++) {
	                p[i6][i5] = p[0][0]+n4*(i5+n5*i6);
			for (i4=0; i4<n4; i4++) {
			        p[i6][i5][i4] = p[0][0][0]+n3*(i4+n4*(i5+n5*i6));
				for (i3=0; i3<n3; i3++) {
				        p[i6][i5][i4][i3] = p[0][0][0][0]
					      +n2*(i3+n3*(i4+n4*(i5+n5*i6)));
					for (i2=0; i2<n2; i2++)
					        p[i6][i5][i4][i3][i2] = 
						      (char*)p[0][0][0][0][0]+
				       size*n1*(i2+n2*(i3+n3*(i4+n4*(i5+n5*i6))));
			        }
		        }
	        }
        }
	return p;
}

/* free a 6-d array */
void free6 (void ******p)
{
        free(p[0][0][0][0][0]);
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, size_t n1)
{
	return (int*)realloc1(v,n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(size_t n1, size_t n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(size_t n1, size_t n2, size_t n3)
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, size_t n1)
{
	return (float*)realloc1(v,n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***)p);
}

/* allocate a 4-d array of floats, added by Zhaobo Meng, 1997 */
float ****alloc4float(size_t n1, size_t n2, size_t n3, size_t n4)
{
        return (float****)alloc4(n1,n2,n3,n4,sizeof(float));
}

/* free a 4-d array of floats, added by Zhaobo Meng, 1997 */
void free4float(float ****p)
{
        free4((void****)p);
}

/* allocate a 5-d array of floats, added by Zhaobo Meng, 1997 */
float *****alloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        return (float*****)alloc5(n1,n2,n3,n4,n5,sizeof(float));
}

/* free a 5-d array of floats, added by Zhaobo Meng, 1997 */
void free5float(float *****p)
{
        free5((void*****)p);
}

/* allocate a 6-d array of floats, added by Zhaobo Meng, 1997 */
float ******alloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6)
{
        return (float******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(float));
}

/* free a 6-d array of floats, added by Zhaobo Meng, 1997 */
void free6float(float ******p)
{
        free6((void******)p);
}

/* allocate a 4-d array of ints, added by Zhaobo Meng, 1997 */
int ****alloc4int(size_t n1, size_t n2, size_t n3, size_t n4)
{
        return (int****)alloc4(n1,n2,n3,n4,sizeof(int));
}

/* free a 4-d array of ints, added by Zhaobo Meng, 1997 */
void free4int(int ****p)
{
        free4((void****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
int *****alloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        return (int*****)alloc5(n1,n2,n3,n4,n5,sizeof(int));
}

/* free a 5-d array of ints, added by Zhaobo Meng, 1997 */
void free5int(int *****p)
{
        free5((void*****)p);
}

/* allocate a 5-d array of chars, added by Zhaobo Meng, 1997 */
unsigned char *****alloc5uchar(size_t n1, size_t n2, size_t n3, 
	size_t n4, size_t n5)
{
        return (unsigned char*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned char));
}

/* free a 5-d array of chars, added by Zhaobo Meng, 1997 */
void free5uchar(unsigned char *****p)
{
        free5((void*****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short *****alloc5ushort(size_t n1, size_t n2, size_t n3,
        size_t n4, size_t n5)
{
        return (unsigned short*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned short));
}

/* free a 5-d array of shorts, added by Zhaobo Meng, 1997 */
void free5ushort(unsigned short *****p)
{
        free5((void*****)p);
}

/* allocate a 6-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short ******alloc6ushort(size_t n1, size_t n2, size_t n3,
        size_t n4, size_t n5, size_t n6)
{
        return (unsigned short******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(unsigned short));
}

/* free a 6-d array of shorts, added by Zhaobo Meng, 1997 */
void free6ushort(unsigned short ******p)
{
        free6((void******)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(size_t n1)
{
	return (double*)alloc1(n1,sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, size_t n1)
{
	return (double*)realloc1(v,n1,sizeof(double));
}


/* free a 1-d array of doubles */
void free1double(double *p)
{
	free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(size_t n1, size_t n2)
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(size_t n1, size_t n2, size_t n3)
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
	free3((void***)p);
}


