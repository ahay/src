#include <rsf.h>
#include "_cjb.h"

/******************************************************************
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
alloc2ushort    allocate a 2-d array of unsigned shorts 
free2ushort     free a 2-d array of unsiged shorts
alloc3ushort    allocate a 3-d array of unsigned shorts 
free3ushort     free a 3-d array of unsiged shorts
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

zero1int        initialize the 1-d int array with zero
zero2int        initialize the 2-d int array with zero
zero3int        initialize the 3-d int array with zero

zero1float      initialize the 1-d float array with zero
zero2float      initialize the 2-d float array with zero
zero3float      initialize the 3-d float array with zero
zero4float      initialize the 4-d float array with zero

zero1double     initialize the 1-d double array with zero
zero2double     initialize the 2-d double array with zero
zero3double     initialize the 3-d double array with zero

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
functions described by Press et al, 1988, Numerical Recipes in C.
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
**************************************************************************/

void *alloc1 (int n1, int size)
/*< allocate a 1-d array >*/
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

void *realloc1(void *v, int n1, int size)
/*< re-allocate a 1-d array >*/
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

void free1 (void *p)
/*< free a 1-d array >*/
{
	free(p);
}

void **alloc2 (int n1, int n2, int size)
/*< allocate a 2-d array >*/
{
	int i2;
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

void free2 (void **p)
/*< free a 2-d array >*/
{
	free(p[0]);
	free(p);
}

void ***alloc3 (int n1, int n2, int n3, int size)
/*< allocate a 3-d array >*/
{
	int i3,i2;
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

void free3 (void ***p)
/*< free a 3-d array >*/
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

void ****alloc4 (int n1, int n2, int n3, int n4, int size)
/*< allocate a 4-d array >*/
{
	int i4,i3,i2;
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

void free4 (void ****p)
/*< free a 4-d array >*/
{
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

void *****alloc5 (int n1, int n2, int n3, int n4, int n5, int size)
/*< The following two functions were added by Zhaobo Meng, Jan. 1997
 allocate a 5-d array >*/
{
	int i5,i4,i3,i2;
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

void free5 (void *****p)
/*< free a 5-d array >*/
{
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

void ******alloc6 (int n1, int n2, int n3, int n4, int n5, int n6,
                  int size)
/*< The following two functions were added by Zhaobo Meng, Jan. 1997
 allocate a 6-d array >*/
{
	int i6,i5,i4,i3,i2;
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

void free6 (void ******p)
/*< free a 6-d array >*/
{
        free(p[0][0][0][0][0]);
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

char *alloc1char(int n1)
/*< allocate a 1-d array of char >*/
{
	return (char*)alloc1(n1,sizeof(char));
}

int *alloc1int(int n1)
/*< allocate a 1-d array of ints >*/
{
	return (int*)alloc1(n1,sizeof(int));
}

int *realloc1int(int *v, int n1)
/*< re-allocate a 1-d array of ints >*/
{
	return (int*)realloc1(v,n1,sizeof(int));
}

void free1char(char *p)
/*< free a 1-d array of char >*/
{
	free1(p);
}

void free1int(int *p)
/*< free a 1-d array of ints >*/
{
	free1(p);
}

char **alloc2char(int n1, int n2)
/*< allocate a 2-d array of char 
 n1: fast dimension; n2: slow dimension >*/
{
	return (char**)alloc2(n1,n2,sizeof(char));
}

int **alloc2int(int n1, int n2)
/*< allocate a 2-d array of ints 
  n1: fast dimension; n2: slow dimension >*/
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

void free2char(char **p)
/*< free a 2-d array of char >*/
{
	free2((void**)p);
}

void free2int(int **p)
/*< free a 2-d array of ints >*/
{
	free2((void**)p);
}

int ***alloc3int(int n1, int n2, int n3)
/*< allocate a 3-d array of ints >*/
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

void free3int(int ***p)
/*< free a 3-d array of ints >*/
{
	free3((void***)p);
}

float *alloc1float(int n1)
/*< allocate a 1-d array of floats >*/
{
	return (float*)alloc1(n1,sizeof(float));
}

float *realloc1float(float *v, int n1)
/*< re-allocate a 1-d array of floats >*/
{
	return (float*)realloc1(v,n1,sizeof(float));
}

void free1float(float *p)
/*< free a 1-d array of floats >*/
{
	free1(p);
}

float **alloc2float(int n1, int n2)
/*< allocate a 2-d array of floats 
  n1: fast dimension; n2: slow dimension >*/
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

void free2float(float **p)
/*< free a 2-d array of floats >*/
{
	free2((void**)p);
}

float ***alloc3float(int n1, int n2, int n3)
/*< allocate a 3-d array of floats >*/
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

void free3float(float ***p)
/*< free a 3-d array of floats >*/
{
	free3((void***)p);
}

float ****alloc4float(int n1, int n2, int n3, int n4)
/*< allocate a 4-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        return (float****)alloc4(n1,n2,n3,n4,sizeof(float));
}

void free4float(float ****p)
/*< free a 4-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        free4((void****)p);
}

float *****alloc5float(int n1, int n2, int n3, int n4, int n5)
/*< allocate a 5-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        return (float*****)alloc5(n1,n2,n3,n4,n5,sizeof(float));
}

void free5float(float *****p)
/*< free a 5-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        free5((void*****)p);
}

float ******alloc6float(int n1, int n2, int n3, int n4, int n5, int n6)
/*< allocate a 6-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        return (float******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(float));
}

void free6float(float ******p)
/*< free a 6-d array of floats, added by Zhaobo Meng, 1997 >*/
{
        free6((void******)p);
}

int ****alloc4int(int n1, int n2, int n3, int n4)
/*< allocate a 4-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        return (int****)alloc4(n1,n2,n3,n4,sizeof(int));
}

void free4int(int ****p)
/*< free a 4-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        free4((void****)p);
}

int *****alloc5int(int n1, int n2, int n3, int n4, int n5)
/*< allocate a 5-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        return (int*****)alloc5(n1,n2,n3,n4,n5,sizeof(int));
}

void free5int(int *****p)
/*< free a 5-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        free5((void*****)p);
}

unsigned char *****alloc5uchar(int n1, int n2, int n3, 
	int n4, int n5)
/*< allocate a 5-d array of chars, added by Zhaobo Meng, 1997 >*/
{
        return (unsigned char*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned char));
}

void free5uchar(unsigned char *****p)
/*< free a 5-d array of chars, added by Zhaobo Meng, 1997 >*/
{
        free5((void*****)p);
}

unsigned short *****alloc5ushort(int n1, int n2, int n3, int n4, int n5)
/*< allocate a 5-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        return (unsigned short*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned short));
}

unsigned short ***alloc3ushort(int n1, int n2, int n3)
/*< allocate a 3-d array of ints, added by Meng, 1997 >*/
{
        return (unsigned short***)alloc3(n1,n2,n3,sizeof(unsigned short));
}

unsigned short **alloc2ushort(int n1, int n2)
/*< allocate a 2-d array of ints, added by Meng, 1997 >*/
{
        return (unsigned short**)alloc2(n1,n2,sizeof(unsigned short));
}

void free5ushort(unsigned short *****p)
/*< free a 5-d array of shorts, added by Zhaobo Meng, 1997 >*/
{
        free5((void*****)p);
}

void free3ushort(unsigned short ***p)
/*< free a 3-d array of shorts, added by Zhaobo Meng, 1997 >*/
{
        free3((void***)p);
}

void free2ushort(unsigned short **p)
/*< free a 2-d array of shorts, added by Zhaobo Meng, 1997 >*/
{
        free2((void**)p);
}

unsigned short ******alloc6ushort(int n1, int n2, int n3,
        int n4, int n5, int n6)
/*< allocate a 6-d array of ints, added by Zhaobo Meng, 1997 >*/
{
        return (unsigned short******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(unsigned short));
}

void free6ushort(unsigned short ******p)
/*< free a 6-d array of shorts, added by Zhaobo Meng, 1997 >*/
{
        free6((void******)p);
}

double *alloc1double(int n1)
/*< allocate a 1-d array of doubles >*/
{
	return (double*)alloc1(n1,sizeof(double));
}

double *realloc1double(double *v, int n1)
/*< re-allocate a 1-d array of doubles >*/
{
	return (double*)realloc1(v,n1,sizeof(double));
}


void free1double(double *p)
/*< free a 1-d array of doubles >*/
{
	free1(p);
}

double **alloc2double(int n1, int n2)
/*< allocate a 2-d array of doubles 
 n1: fast dimension; n2: slow dimension >*/
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

void free2double(double **p)
/*< free a 2-d array of doubles >*/
{
	free2((void**)p);
}

double ***alloc3double(int n1, int n2, int n3)
/*< allocate a 3-d array of doubles >*/
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

void free3double(double ***p)
/*< free a 3-d array of doubles >*/
{
	free3((void***)p);
}
