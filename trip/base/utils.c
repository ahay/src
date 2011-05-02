/* Utilities. */
/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* 
utils.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "mm_malloc.h"

#include "_cstd.h"
#include "_usempi.h"
/*^*/

#ifndef _sf_utils_h

/*----------------------------------------------------------------------------*/
/* inline keyword. */
#ifndef __cplusplus
    #if __STDC_VERSION__ >= 199901L
    #else
        #define inline
    #endif
#endif
/* restrict keyword. */
#if __STDC_VERSION__ >= 199901L
#else
    #define restrict
#endif
/*----------------------------------------------------------------------------*/
/* 
Data types.
*/
#define DT_CSTRING     1
#define DT_CHAR        2
#define DT_SHORT       3
#define DT_INT         4
#define DT_LONG        5
#define DT_FLOAT       6
#define DT_DOUBLE      7
#define DT_USHORT      9
#define DT_UINT       10
#define DT_ULONG      11
/*----------------------------------------------------------------------------*/
/**
Real type and consts.- legal values are DT_DOUBLE and DT_FLOAT.
*/
/*#define DT_REAL DT_DOUBLE*/  
#define DT_REAL DT_FLOAT
/*----------------------------------------------------------------------------*/
/* DO NOT CHANGE. */
#if   DT_REAL == DT_DOUBLE
/** real type */
    typedef double ireal;
    #define REAL_NAN 0.0         /* TODO: better than 0 */
    #define REAL_ZERO 0.0
    #define REAL_ONE 1.0
    #define REAL_EPS DBL_EPSILON
    #define REAL_MIN DBL_MIN
    #define REAL_MAX DBL_MAX
	#define IWAVE_MPI_REAL MPI_DOUBLE
#elif DT_REAL == DT_FLOAT
/** real type */
    typedef float ireal;
    #define REAL_NAN 0.0f        /* TODO: better than 0 */
    #define REAL_ZERO 0.0f
    #define REAL_ONE 1.0f
    #define REAL_EPS FLT_EPSILON
    #define REAL_MIN FLT_MIN
    #define REAL_MAX FLT_MAX
	#define IWAVE_MPI_REAL MPI_FLOAT
#else
	#error REAL TYPE UNDEFINED.
#endif
/*^*/

/*----------------------------------------------------------------------------*/
/**
Maximum number of arrays in any domain.
*/
#define RDOM_MAX_NARR 20
/*----------------------------------------------------------------------------*/
/**
Maximum number of dimensions in any array. If set > 3, note
that RARR_MAX_3NDIM (below) will require redefinition.
*/
#define RARR_MAX_NDIM 3
/*^*/

/******************* begin formerly in iwave.h - moved 26.01.11 WWS ***********/
/*
  Number of dimensions. 1D/2D/3D only. 
  CHANGE BEFORE COMPILING. 
  IWAVE_NDIM SHOULD BE >= RARR_MAX_NDIM.
*/
#define IWAVE_NDIM 3   /* 1 | 2 | 3 */
/*^*/
  
#if IWAVE_NDIM > RARR_MAX_NDIM
#error NOT ENOUGH DIMENSIONS IN THE RARRAY: WAVE_NDIM > RARR_MAX_NDIM.
#endif
/*----------------------------------------------------------------------------*/
/*
  3^IWAVE_NDIM and number of neighbors for 1D/2D/3D.
  Note: lranks index is computed via ex_ind2n(..), 
  lranks[IWAVE_NNEI] contains this processor's rank.
*/
#if   IWAVE_NDIM == 1
#define IWAVE_3NDIM 3
#elif IWAVE_NDIM == 2
#define IWAVE_3NDIM 9
#elif IWAVE_NDIM == 3
#define IWAVE_3NDIM 27
#else
#error IWAVE_3NDIM UNDEFINED.
#endif
/*^*/

#define IWAVE_NNEI (IWAVE_3NDIM - 1)
/*^*/

/********************* end formerly in iwave.h - moved 26.01.11 WWS ***********/

/*----------------------------------------------------------------------------*/
/*
Indices of coordinate axes. Permutation of {0,1,2}.
TO CHANGE BY USER.
*/
#define IX 0
#define IY 1
#define IZ 2
/*----------------------------------------------------------------------------*/
/*
Flag (0/1) to dump pointers when dumping an array.
TO CHANGE BY USER.
*/
#define RARR_DUMP_POINTERS 0
/*----------------------------------------------------------------------------*/
/**
\section IPNT
Integer arrays of space dimension - index into grid
*/
typedef int IPNT[RARR_MAX_NDIM];
/**
\section RPNT
Real arrays of space dimension - grid coordinates, vector field values.
*/
typedef ireal RPNT[RARR_MAX_NDIM];
/*^*/

#if   RARR_MAX_NDIM == 1
    static const IPNT IPNT_0 = {0};
    static const IPNT IPNT_1 = {1};
    static const RPNT RPNT_0 = {REAL_ZERO};
	static const RPNT RPNT_1 = {REAL_ONE};
    #define RARR_MAX_3NDIM     3
#elif RARR_MAX_NDIM == 2
    static const IPNT IPNT_0 = {0, 0};
    static const IPNT IPNT_1 = {1, 1};
    static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO};
	static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE};
    #define RARR_MAX_3NDIM     9
#elif RARR_MAX_NDIM == 3
    static const IPNT IPNT_0 = {0, 0, 0};
    static const IPNT IPNT_1 = {1, 1, 1};
    static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO, REAL_ZERO};
	static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE, REAL_ONE};
    #define RARR_MAX_3NDIM     27
#else
	#error IPNT/RPNT CONSTANTS UNDEFINED.
#endif
/*^*/

/* compatibility */
#define _IPNT IPNT
#define _RPNT RPNT
/*----------------------------------------------------------------------------*/
/*
Flag to perform bounds checks.
TO CHANGE BY USER.
*/
#define CHECK_BOUNDS
/*#undef CHECK_BOUNDS*/
/*----------------------------------------------------------------------------*/
/**
String with size. Has size n (0,...,n-1) and extra terminating null character 
at n-th position. May have null characters in the middle.
*
long n  :  string length (not including null terminator).
char *s :  string pointer (has additional null terminator).
*/
typedef struct
{
    long n;
    char *s;
} SIZEDSTRING;
/*----------------------------------------------------------------------------*/
/*
Separator and quote symbols for parser.
TO CHANGE BY USER.
*/
#define PS_SEP '='     /* separator symbol */
#define PS_QUO '"'     /* quote symbol */
/*^*/

#define iwave_min(a, b) ((a) < (b) ? (a) : (b))
#define iwave_max(a, b) ((a) > (b) ? (a) : (b))
/*^*/

/*@{*/
#define E_SUCCESS          0     /**< no error */
#define E_INTERNAL      1000     /**< internal error, should not happen */
#define E_OTHER         2000     /**< some other error */
#define E_ALLOC            1     /**< cannot (re/m/c)alloc, etc */
#define E_BADINPUT         2     /**< incorrect input data */
#define E_OUTOFBOUNDS      3     /**< if array bound check failed */
#define E_BADINDEX         4     /**< wrong index */
#define E_BADARRINDEX      5     /**< wrong array index in the domain */
#define E_BADDIMINDEX      6     /**< wrong dimension index in the array */
#define E_FILE             7     /**< problem with file operations */
#define E_FILEOPEN         8     /**< problem opening file */
#define E_MPI              9     /**< mpi function returned error code */
#define E_DOMAINDECOMP    10     /**< could not do domain decomposition */
#define E_PARSE           11     /**< problem with parsing parameter file */
#define E_PARSENONAME     12     /**< name not found */
#define E_PARSENOVALUE    13     /**< value not found */
#define E_PARSECONVERT    14     /**< cannot convert */
#define E_ALREADYALLOC    15     /**< allocating already allocated array */
#define E_RANGE           -2     /**< out of range */
#define E_OVERFLOW        -3     /**< out of range - overflow */
#define E_UNDERFLOW       -4     /**< out of range - underflow */
#define E_NOTIMESTEP    -100     /**< no time step needed for array */
#define E_NOTINGRID     -200     /**< process with these indices not in domain decomposition grid */ 
/*@}*/
/*^*/

#define SEAMX_BIG_ENDIAN    0
#define SEAMX_LITTLE_ENDIAN 1
/*^*/

#endif

/*----------------------------------------------------------------------------*/

static int m_rk = 0;
static int m_sz = 1;
static MPI_Comm m_cm = 0;
static FILE *m_stream = 0;

/* added 31.01.09 WWS
 */
static int m_grk = 0;
static int m_gsz = 1;
static MPI_Comm m_gcm = 0;
static int m_ts = 1;

/* added 21.11.10 WWS
 */
static int m_g = 0;
static int m_ng = 1;

/* added 17.01.11 WWS
 */
static MPI_Comm m_rcm = 0;

/*----------------------------------------------------------------------------*/

void* usermalloc_(size_t size)
/*<    return malloc(size); >*/
{
    return _mm_malloc(size, 16);
}
/*----------------------------------------------------------------------------*/

void userfree_(void *ptr)
/*<    free(ptr); >*/
{
    _mm_free(ptr);
}
/*----------------------------------------------------------------------------*/

int* IASN(IPNT l, const IPNT r)
/*< ** Assignment of IPNT r to IPNT l >*/
{
    if ( l == r ) return l;
    return (int*)memcpy((void*)l, (const void*)r, RARR_MAX_NDIM * sizeof(int));
}
/*----------------------------------------------------------------------------*/

ireal* RASN(RPNT l, const RPNT r)
/*< ** Assignment of RPNT r to RPNT l >*/
{
    if ( l == r ) return l;
    return (ireal*)memcpy((void*)l, (const void*)r, RARR_MAX_NDIM * sizeof(ireal));
}
/*----------------------------------------------------------------------------*/

int gen_3n1(int ndim, int *n)
/*< Numbering functions.
Convert {-1,0,1}^ndim to a linear index and back. 
Note: {0,0,0} always gets last linear index (3^ndim - 1).

gen_3n1   :  number of linear indices - 1 (IMPORTANT:  returns n = 3^ndim - 1).
gen_n2pnt :  linear index to cartesian.
gen_pnt2n :  cartesian index to linear.
  
int ndim   :  number of dimensions.
IPNT p     :  cartesian index, each p[i] is from {-1,0,1} set.
int i      :  linear index from 0 to 3^ndim - 1 (including).
int return :  error code.
>*/
{
    int data[3] = {2, 8, 26};
    
    if ( (unsigned int)(--ndim) >= (unsigned int)3 ) return E_BADINPUT;
    
    *n = data[ndim];
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int gen_i2pnt(int ndim, int i, IPNT p)
/*< Numbering functions.
Convert {-1,0,1}^ndim to a linear index and back. 
Note: {0,0,0} always gets last linear index (3^ndim - 1).

gen_3n1   :  number of linear indices - 1 (IMPORTANT:  returns n = 3^ndim - 1).
gen_n2pnt :  linear index to cartesian.
gen_pnt2n :  cartesian index to linear.
  
int ndim   :  number of dimensions.
IPNT p     :  cartesian index, each p[i] is from {-1,0,1} set.
int i      :  linear index from 0 to 3^ndim - 1 (including).
int return :  error code.
>*/
{
    int idim, k, imax;
    
    if ( gen_3n1(ndim, &imax) ) return E_BADINPUT;
    if ( (unsigned int)i > (unsigned int)imax ) return E_BADINPUT;
    
    for ( idim = 0; idim < ndim; ++idim )
    {
        k = i % 3;
        i /= 3;
        if      ( k == 0 ) p[idim] = -1;
        else if ( k == 2 ) p[idim] = 0;
        else if ( k == 1 ) p[idim] = 1;
        else return E_INTERNAL;
    }

    return 0;
}
/*----------------------------------------------------------------------------*/

int gen_pnt2i(int ndim, const IPNT p, int *i)
/*< Numbering functions.
Convert {-1,0,1}^ndim to a linear index and back. 
Note: {0,0,0} always gets last linear index (3^ndim - 1).

gen_3n1   :  number of linear indices - 1 (IMPORTANT:  returns n = 3^ndim - 1).
gen_n2pnt :  linear index to cartesian.
gen_pnt2n :  cartesian index to linear.
  
int ndim   :  number of dimensions.
IPNT p     :  cartesian index, each p[i] is from {-1,0,1} set.
int i      :  linear index from 0 to 3^ndim - 1 (including).
int return :  error code.
>*/
{
    int idim, k;
    
    if ( (unsigned int)(--ndim) >= (unsigned int)3 ) return E_BADINPUT;

    *i = 0;
    for ( idim = ndim; idim >= 0; --idim )
    {
        if      ( p[idim] == -1 ) k = 0;
        else if ( p[idim] ==  0 ) k = 2; /* code 0 as 2, so [0,...,0] is last */
        else if ( p[idim] ==  1 ) k = 1;
        else return E_BADINPUT;
        *i = 3 * (*i) + k;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

void storeRank(int rk)
/*< store rank >*/
{
	m_rk = rk;
}
/*----------------------------------------------------------------------------*/

int retrieveRank()
/*< retrieve rank >*/
{
	return m_rk;
}
/*----------------------------------------------------------------------------*/

void storeSize(int sz)
/*< store size >*/
{
	m_sz = sz;
}
/*----------------------------------------------------------------------------*/

int retrieveSize()
/*< retrieve size >*/
{
	return m_sz;
}
/*----------------------------------------------------------------------------*/

void storeComm(MPI_Comm cm)
/*< store communicator >*/
{
	m_cm = cm;
}
/*----------------------------------------------------------------------------*/

MPI_Comm retrieveComm()
/*< retrieve communicator >*/
{
	return m_cm;
}
/*----------------------------------------------------------------------------*/

void storeOutstream(FILE *stream)
/*< store output stream >*/
{
	m_stream = stream;
}
/*----------------------------------------------------------------------------*/

FILE* retrieveOutstream()
/*< retrieve output stream >*/
{
	return m_stream;
}
/* added 31.01.09 WWS */
/*----------------------------------------------------------------------------*/

void storeGlobalRank(int rk)
/*< store global rank >*/
{
	m_grk = rk;
}
/*----------------------------------------------------------------------------*/

int retrieveGlobalRank()
/*< retrieve global rank >*/
{
	return m_grk;
}
/*----------------------------------------------------------------------------*/

void storeGlobalSize(int sz)
/*< store global size >*/
{
	m_gsz = sz;
}
/*----------------------------------------------------------------------------*/

int retrieveGlobalSize()
/*< retrieve global size >*/
{
	return m_gsz;
}
/*----------------------------------------------------------------------------*/

void storeGlobalComm(MPI_Comm cm)
/*< store global communicator >*/
{
	m_gcm = cm;
}
/*----------------------------------------------------------------------------*/

MPI_Comm retrieveGlobalComm()
/*< retrieve global communicator >*/
{
	return m_gcm;
}
/*----------------------------------------------------------------------------*/

void storeThreadSupp(int ts)
/*< store thread >*/
{
	m_ts = ts;
}
/*----------------------------------------------------------------------------*/

int retrieveThreadSupp()
/*< retrieve thread >*/
{
	return m_ts;
}

/* added 21.11.10 WWS */
/*----------------------------------------------------------------------------*/
void storeGroupID(int g) 
/*< >*/
{ m_g=g; }

/*----------------------------------------------------------------------------*/
int retrieveGroupID() 
/*< >*/
{ return m_g; }

/*----------------------------------------------------------------------------*/
void storeNumGroups(int ng) 
/*< >*/
{ m_ng = ng; }

/*----------------------------------------------------------------------------*/
int retrieveNumGroups() 
/*< >*/
{ return m_ng; }

/* added 17.01.11 WWS */
/*----------------------------------------------------------------------------*/
void storeRemComm(MPI_Comm cm) 
/*< >*/
{ m_rcm = cm; }

/*----------------------------------------------------------------------------*/
MPI_Comm retrieveRemComm() 
/*< >*/
{ return m_rcm; }

/*----------------------------------------------------------------------------*/
int getMachineEndianness()
/*< get endianness >*/
{
	short int word = 0x0001;
	char *byte = (char *) &word;
	return (byte[0] ? SEAMX_LITTLE_ENDIAN : SEAMX_BIG_ENDIAN);
}
/*----------------------------------------------------------------------------*/

void swapBytes(unsigned char *arr, int arrsize, int atomsize)
/*< get endianness >*/
{
	register int i, j, k, atomsize1;
	register unsigned char tmp;
	
	atomsize1 = atomsize - 1;
	
	for ( i = 0; i < arrsize; ++i )
	{
		j = 0;
		k = atomsize1;
		while (j < k)
		{
			tmp = arr[j];
			arr[j++] = arr[k];
			arr[k--] = tmp;
		}
		arr += atomsize;
	}
}
/*----------------------------------------------------------------------------*/

