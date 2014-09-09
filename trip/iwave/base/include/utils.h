/** \file
    Author: Igor Terentyev.
    ********************************************************************************
    Utilities:

    Constants.
    Data type [ireal].
    Data type flags.
    Point arrays.
    Aligned malloc.
    Separator, delimiter for parser.
    All error codes.
*/
/*============================================================================*/

#ifndef __UTILS_H_
#define __UTILS_H_
/*----------------------------------------------------------------------------*/

#include "cstd.h"
#include "usempi.h"

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
#define ireal double
#define REAL_NAN 0.0         /* TODO: better than 0 */
#define REAL_ZERO 0.0
#define REAL_ONE 1.0
#define REAL_EPS DBL_EPSILON
#define REAL_MIN DBL_MIN
#define REAL_MAX DBL_MAX
#define IWAVE_MPI_REAL MPI_DOUBLE
#elif DT_REAL == DT_FLOAT
#define ireal float
/**< ** real type >*/
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

/*----------------------------------------------------------------------------*/
/**
   Maximum number of arrays in any domain.
*/
#define RDOM_MAX_NARR 30
/*----------------------------------------------------------------------------*/
/**
   Maximum number of dimensions in any array. If set > 3, note
   that RARR_MAX_3NDIM (below) will require redefinition.
*/
#define RARR_MAX_NDIM 7

/******************* begin formerly in iwave.h - moved 26.01.11 WWS ***********/
/*
  Number of SPATIAL dimensions. 1D/2D/3D only. 
  CHANGE BEFORE COMPILING. 
  IWAVE_NDIM SHOULD BE >= RARR_MAX_NDIM.
*/
#define IWAVE_NDIM 3   /* 1 | 2 | 3 */
  
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

#define IWAVE_NNEI (IWAVE_3NDIM - 1)

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
   Integer arrays of at least space dimension - index into grid
*/
#define IARR_MAX_NDIM 7

#if RARR_MAX_NDIM > IARR_MAX_NDIM
#error NOT ENOUGH DIMENSIONS IN THE IPNT: RARR_MAX_NDIM > IARR_MAX_NDIM.
#endif
typedef int IPNT[IARR_MAX_NDIM];
/**
   \section RPNT
   Real arrays of space dimension - grid coordinates, vector field values.
*/
typedef ireal RPNT[RARR_MAX_NDIM];

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
#elif RARR_MAX_NDIM == 4
static const IPNT IPNT_0 = {0, 0, 0, 0};
static const IPNT IPNT_1 = {1, 1, 1, 1};
static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO};
static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE};
#define RARR_MAX_3NDIM     81
#elif RARR_MAX_NDIM == 5
static const IPNT IPNT_0 = {0, 0, 0, 0, 0};
static const IPNT IPNT_1 = {1, 1, 1, 1, 1};
static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO};
static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE};
#define RARR_MAX_3NDIM     243
#elif RARR_MAX_NDIM == 6
static const IPNT IPNT_0 = {0, 0, 0, 0, 0, 0};
static const IPNT IPNT_1 = {1, 1, 1, 1, 1, 1};
static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO};
static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE};
#define RARR_MAX_3NDIM     729
#elif RARR_MAX_NDIM == 7
static const IPNT IPNT_0 = {0, 0, 0, 0, 0, 0, 0};
static const IPNT IPNT_1 = {1, 1, 1, 1, 1, 1, 1};
static const RPNT RPNT_0 = {REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO, REAL_ZERO};
static const RPNT RPNT_1 = {REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE, REAL_ONE};
#define RARR_MAX_3NDIM     2187
#else
#error IPNT/RPNT CONSTANTS UNDEFINED.
#endif



/** Assignment of IPNT r to IPNT l */
int* IASN(IPNT l, const IPNT r);
/** Assignment of RPNT r to RPNT l */
ireal* RASN(RPNT l, const RPNT r);

/* compatibility */
#define _IPNT IPNT
#define _RPNT RPNT
/*----------------------------------------------------------------------------*/
/*
  Flag to perform bounds checks.
  TO CHANGE BY USER.
*/
#define CHECK_BOUNDS
/* #undef CHECK_BOUNDS */
/*----------------------------------------------------------------------------*/
/**
   String with size. Has size n (0,...,n-1) and extra terminating null character 
   at n-th position. May have null characters in the middle.

   long n  :  string length (not including null terminator).
   char *s :  string pointer (has additional null terminator).
*/
typedef struct s_SIZEDSTRING
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
/*----------------------------------------------------------------------------*/
void* usermalloc_(size_t size);
/**<
   Modified malloc/free functions.
   Currently implemented as 16 byte aligned pointers.
   TO CHANGE BY USER in utils.c.

   size_t size       :  size in bytes.
   void * return/ptr :  pointer to the allocated memory.
   >*/
void userfree_(void *ptr);
/*----------------------------------------------------------------------------*/
/*
  Min, max macros.
*/

#define iwave_min(a, b) ((a) < (b) ? (a) : (b))
#define iwave_max(a, b) ((a) > (b) ? (a) : (b))
#define iwave_abs(a) ((a) < REAL_ZERO ? (-a) : (a))

/*----------------------------------------------------------------------------*/
/** \defgroup error Error codes
    Error codes used throughout IWAVE.
*/
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

/*----------------------------------------------------------------------------*/
/*
  Numbering functions.
  Convert {-1,0,1}^ndim to a linear index and back. 
  Note: {0,0,0} always gets last linear index (3^ndim - 1).

  gen_3n1   :  number of linear indices - 1 (IMPORTANT:  returns n = 3^ndim - 1).
  gen_n2pnt :  linear index to cartesian.
  gen_pnt2n :  cartesian index to linear.
  
  int ndim   :  number of dimensions.
  IPNT p     :  cartesian index, each p[i] is from {-1,0,1} set.
  int i      :  linear index from 0 to 3^ndim - 1 (including).
  int return :  error code.
*/
int gen_3n1(int ndim, int *n);
int gen_i2pnt(int ndim, int i, IPNT p);
int gen_pnt2i(int ndim, const IPNT p, int *i);
/*----------------------------------------------------------------------------*/

void storeRank(int rk);
int retrieveRank();

void storeSize(int sz);
int retrieveSize();

void storeComm(MPI_Comm cm);
MPI_Comm retrieveComm();

void storeOutstream(FILE *stream);
FILE* retrieveOutstream();

/*
 * added 31.01.09 WWS
 */
void storeGlobalRank(int rk);
int retrieveGlobalRank();

void storeGlobalSize(int sz);
int retrieveGlobalSize();

void storeGlobalComm(MPI_Comm cm);
MPI_Comm retrieveGlobalComm();

void storeThreadSupp(int ts);
int retrieveThreadSupp();

/*
 * added 21.11.10 WWS
 */
void storeGroupID(int);
int retrieveGroupID();

void storeNumGroups(int);
int retrieveNumGroups();

/*
 * added 17.01.11 WWS
 */
void storeRemComm(MPI_Comm cm);
MPI_Comm retrieveRemComm();

/*----------------------------------------------------------------------------*/
/* 
   Determine machine endiannes runtime
*/
#define SEAMX_BIG_ENDIAN    0
#define SEAMX_LITTLE_ENDIAN 1

int getMachineEndianness();

/*
  Swap bytes for each number in the array.
  Array has [arrsize] elements, each [atomsize] bytes long.
*/
void swapBytes(unsigned char *arr, int arrsize, int atomsize);

#endif
