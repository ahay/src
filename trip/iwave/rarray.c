/* Dimension information, array type and its operations. */
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
rarray.c
Igor Terentyev.
*/
/*============================================================================*/

#include <trip/base.h>
/*^*/

#include "rarray.h"

#include "exchangeinfo.h"
/*^*/

#ifndef _sf_rarray_h

typedef struct
{
  /** the working (computational virtual) array size */
  int n;    
  /** global start index of the working (computational virtual) array */
  int gs;   
  /** global end index of the working (computational virtual) array */
  int ge;   
  /** size of the allocated array, [used for stride] */
  int n0;   
  /** global start index of the allocated array [used for bounds checks] */
  int gs0;  
  /** global end index of the allocated array [used for bounds checks] */
  int ge0;  
} INFODIM;
/*^*/

typedef struct
{
  /** number of dimensions, e.g.\ 1,2,3 */
  int ndim; 
  /** beginning pointer of the working (computational virtual) array */
  ireal *_s;
  /** beginning pointer of the allocated array */            
  ireal *_s0;      
  /** 
   * dimension information. RARR_MAX_NDIM: 
   * maximum dimension number 
   */
  INFODIM _dims[RARR_MAX_NDIM];  
} RARR;
/*^*/

typedef int (RA_CREATE_FUN)(RARR *arr, int ndim, IPNT v1, IPNT v2);
/*^*/

typedef int (RA_SET_FUN)(RARR *arr, const IPNT v1, const IPNT v2);
/*^*/

#endif

/*----------------------------------------------------------------------------*/
/*
Out of bounds message.
*/
static void ra_oob(const RARR *arr, int idim, int ind, const char *funname)
{
    fprintf(stderr, "%s OUT OF BOUNDS: dim = %d, ind = %d\n", funname, idim, ind);
    ra_dump(arr, stderr);
}
/*----------------------------------------------------------------------------*/

int ra_setnull(RARR *arr)
/*< Set all fields in an array to zeros (NO DEALLOCATION). >*/
{
    memset((void*)arr, 0, sizeof(RARR));
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_create(RARR *arr, int ndim, IPNT gs, IPNT ge)
/*<* 
 * Create array (given gs and n). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0
>*/
{
    int err;
    
    err = ra_declare(arr, ndim, gs, ge);
    if ( err ) return err;
    
    err = ra_allocate(arr);
    if ( err )
    {    
        ra_setnull(arr);
        return err;
    }

    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_create_s(RARR *arr, int ndim, IPNT gs, IPNT n)
/*<* 
 * Create array (given gs and n). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0
>*/
{
    IPNT ge; 
    int d;
    
    ndim = iwave_min(ndim, RARR_MAX_NDIM);
    for ( d = 0; d < ndim; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_create(arr, ndim, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_create_e(RARR *arr, int ndim, IPNT ge, IPNT n)
/*<* 
 * Create array (given ge and n). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0
>*/
{
    IPNT gs;
    int d;
    
    ndim = iwave_min(ndim, RARR_MAX_NDIM);
    for ( d = 0; d < ndim; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_create(arr, ndim, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare_s(RARR *arr, int ndim, IPNT gs, IPNT n)
/*<* 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0
>*/
{
    IPNT ge; 
    int d;
    
    ndim = iwave_min(ndim, RARR_MAX_NDIM);
    for ( d = 0; d < ndim; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_declare(arr, ndim, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare_e(RARR *arr, int ndim, IPNT ge, IPNT n)
/*<* 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0
>*/
{
    IPNT gs;
    int d;
    
    ndim = iwave_min(ndim, RARR_MAX_NDIM);
    for ( d = 0; d < ndim; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_declare(arr, ndim, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare(RARR *arr, int ndim, IPNT gs, IPNT ge)
/*<* 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @return 0
>*/
{
    int d;
    INFODIM *dim;                      /* pointer to current dimension */
    
    ra_setnull(arr);                   /* empty array */
    if ( (unsigned int)ndim > (unsigned int)RARR_MAX_NDIM ) return E_BADINPUT;

    for ( d = 0; d < ndim; ++d )
    {
        dim = arr->_dims + d;
        dim->gs = dim->gs0 = gs[d];
        dim->ge = dim->ge0 = ge[d];
        dim->n  = dim->n0  = dim->ge - dim->gs + 1;

        if ( dim->n0 < 0 )             /* check negative size */
        {
            ra_setnull(arr);
            return E_BADINPUT;
        }
    }
    
    arr->ndim = ndim;

    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_allocate(RARR *arr)
/*< Allocate array. >*/
{
    int d;
    long size;                          /* size to allocate */
    
    if ( arr->_s0 != NULL ) return E_ALREADYALLOC;
    
    size = 1L;                          /* allocation size */
    for ( d = 0; d < arr->ndim; ++d ) size *= (long)(arr->_dims[d].n0);
    
    if ( size > 0L )                    /* allocate memory */
    {
        arr->_s0 = (ireal*)usermalloc_(size * sizeof(ireal));
        if ( arr->_s0 == NULL ) return E_ALLOC;
    }
    
    arr->_s = arr->_s0;

    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_destroy(RARR *arr)
/*< Destroy array (STORAGE DEALLOCATION). >*/ 
{
  /*
  fprintf(stderr,"ra_destroy: arr->_s0 = %d\n",arr->_s0);
  */
  userfree_(arr->_s0);
  ra_setnull(arr);
  /*
  fprintf(stderr,"ra_destroy exit: arr->_s0 = %d\n",arr->_s0);
  */
  return 0;
}
/*----------------------------------------------------------------------------*/

int ra_offset(RARR *arr, const IPNT os, const IPNT oe)
/*<*
 * Reset the working (computational virtual) array (given oe and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  os  - (IPNT) start index offsets (forward) of the working (computational virtual) array in all dimensions.
 * @param[in]  oe  - (IPNT) end index offsets (backward) of the working (computational virtual) array in all dimensions.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 >*/
{
    IPNT gs, ge;
    int d;
    
    for ( d = 0; d < arr->ndim; ++d )
    {
        gs[d] = arr->_dims[d].gs + os[d];
        ge[d] = arr->_dims[d].ge - oe[d];
    }
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_offset_s(RARR *arr, const IPNT os, const IPNT n)
/*<*
 * Reset the working (computational virtual) array (given os and n) (NO STORAGE ALLOCATION).
 *
 * new_gs = gs + os, new_ge = ge - oe.
 * @param[out] arr - (RARR *)
 * @param[in]  os  - (IPNT) start index offsets (forward) of the working (computational virtual) array in all dimensions.
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT gs, ge;
    int d;
    
    for ( d = 0; d < arr->ndim; ++d )
    {
        gs[d] = arr->_dims[d].gs + os[d];
        ge[d] = gs[d] + n[d] - 1;
    }
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_offset_e(RARR *arr, const IPNT oe, const IPNT n)
/*<*
 * Reset the working (computational virtual) array (given oe and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  oe  - (IPNT) end index offsets (backward) of the working (computational virtual) array in all dimensions.
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT gs, ge;
    int d;
    
    for ( d = 0; d < arr->ndim; ++d )
    {
        ge[d] = arr->_dims[d].ge - oe[d];
        gs[d] = ge[d] - n[d] + 1;
    }
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_greset(RARR *arr, const IPNT gs, const IPNT ge)
/*<*
 * Reset the working (computational virtual) array (given gs and ge) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  gs  - (IPNT) global start indices of the working (computational virtual) array in all dimensions.
 * @param[in]  ge  - (IPNT) global end indices of the working (computational virtual) array in all dimensions.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    int d;
    INFODIM *dim;                      /* pointer to current dimension */
    long pshift = 0L;                  /* pointer shift */
    RARR arrold = *arr;                /* remember old to restore if error */
    
    /* cycle backwards to compute pshift correctly */
    for ( d = arr->ndim - 1; d >= 0; --d ) 
    {
        dim = arr->_dims + d;
        dim->gs = gs[d];
        dim->ge = ge[d];
        dim->n = dim->ge - dim->gs + 1;
        if ( (dim->n < 0) || (dim->gs < dim->gs0) || (dim->ge > dim->ge0) )
        {
	  fprintf(stderr,"n=%d gs=%d gs0=%d ge=%d ge0=%d\n",dim->n,dim->gs,dim->gs0,dim->ge,dim->ge0);
            *arr = arrold;
            return E_OUTOFBOUNDS;
        }
        pshift = pshift * (long)(dim->n0) + (long)(dim->gs - dim->gs0);
    }
    
    if ( arr->_s != NULL ) arr->_s = arr->_s0 + pshift;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_greset_s(RARR *arr, const IPNT gs, const IPNT n)
/*<*
 * Reset the working (computational virtual) array (given gs and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  gs  - (IPNT) global start indices of the working (computational virtual) array in all dimensions
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT ge;
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_greset_e(RARR *arr, const IPNT ge, const IPNT n)
/*<*
 * Reset the working (computational virtual) array (given ge and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  ge  - (IPNT) global end indices of the working (computational virtual) array in all dimensions
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT gs; 
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_dump(const RARR *arr, FILE* stream)
/*< Dump array information. >*/
{
    int d;
    const INFODIM *dim;                      /* pointer to current dimension */
    
    fprintf(stream, "array dump: ndim = %d\n", arr->ndim);

    for ( d = 0; d < arr->ndim; ++d ) 
    {
        dim = arr->_dims + d;
        fprintf( stream, "dim %d:  gs = %d  ge = %d  n = %d  gs0 = %d  ge0 = %d  n0 = %d\n", 
                 d, dim->gs, dim->ge, dim->n, dim->gs0, dim->ge0, dim->n0 );
    }    
    if ( RARR_DUMP_POINTERS )
        fprintf(stream, "_s0 = %p  _s = %p\n", (void*)(arr->_s0), (void*)(arr->_s));
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_fprintslice(RARR *arr, const char *path, int idim, int li)
/*<*
 * Output a slice of the working (computational virtual) array to a stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "w");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = ra_printslice(arr, stream, idim, li);
    fclose(stream);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ra_printslice(RARR *arr, FILE* stream, int idim, int li)
/*<*
 * Output a slice of the working (computational virtual) array to a stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT _li, n;
    int err, d1, d2;
    
    err = ra_checkbound(arr, idim, li);
    if ( err ) return err;
        
    ra_size(arr, n);
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            fprintf(stream, "%+11.3e\n", ra_get(arr, &li));
            break;
        
        case 2:
            _li[idim] = li;
            d1 = 1 - idim;
            for ( _li[d1] = 0; _li[d1] < n[d1]; ++_li[d1] ) 
                fprintf(stream, "%+11.3e ", ra_get(arr, _li));
            fprintf(stream, "\n");    
            break;
            
        case 3:
            _li[idim] = li;
            d1 = (idim > 0) ? 0 : 1;
            d2 = (idim < 2) ? 2 : 1;
            for ( _li[d2] = 0; _li[d2] < n[d2]; ++_li[d2] )
            {
                for ( _li[d1] = 0; _li[d1] < n[d1]; ++_li[d1] ) 
                    fprintf(stream, "%+11.3e ", ra_get(arr, _li));
                fprintf(stream, "\n");    
            }
            break;
        
        default:
            fprintf(stream, "output of arrays of more than 3 dimensions not supported.\n");
            return E_OTHER;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_fwriteslice(RARR *arr, const char *path, int idim, int li)
/*<*
 * Output a slice of the working (computational virtual) array to a file.
 *
 * Format: binary.
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "wb");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = ra_writeslice(arr, stream, idim, li);
    fclose(stream);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ra_writeslice(RARR *arr, FILE* stream, int idim, int li)
/*<*
 * Output a slice of the working (computational virtual) array to a file.
 *
 * Format: binary.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT _li, n;
    int err;
    ireal r;
    
    err = ra_checkbound(arr, idim, li);
    if ( err ) return err;
        
    ra_size(arr, n);
    _li[idim] = li;
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            r = ra_get(arr, _li);
            fwrite(&r, sizeof(ireal), 1, stream);
            break;
        #if RARR_MAX_NDIM > 1
        case 2:
            if ( idim != 0 )
                fwrite(arr->_s + li * arr->_dims[0].n0, sizeof(ireal), n[0], stream);
            else
                for ( _li[1] = 0; _li[1] < n[1]; ++_li[1] ) 
                {
                    r = ra_get(arr, _li);
                    fwrite(&r, sizeof(ireal), 1, stream);
                }
            break;
        #endif
        #if RARR_MAX_NDIM > 2     
        case 3:
            if ( idim != 0 )
                for ( _li[3 - idim] = 0; _li[3 - idim] < n[3 - idim]; ++_li[3 - idim] ) 
                    fwrite(arr->_s + (_li[1] + _li[2] * arr->_dims[1].n0) * arr->_dims[0].n0, 
                           sizeof(ireal), n[0], stream);
            else
                for ( _li[2] = 0; _li[2] < n[2]; ++_li[2] )
                    for ( _li[1] = 0; _li[1] < n[1]; ++_li[1] ) 
                    {
                        r = ra_get(arr, _li);
                        fwrite(&r, sizeof(ireal), 1, stream);
                    }
            break;
        #endif
        default:
            return E_OTHER;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_fprint(RARR *arr, const char *path)
/*<*
 * Output the working (computational virtual) array to a file
 *
 * Format: formatted ASCII
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "w");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = ra_print(arr, stream);
    fclose(stream);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ra_print(RARR *arr, FILE* stream)
/*<*
 * Output the working (computational virtual) array to stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT li, n;
    
    ra_size(arr, n);
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            for ( li[0] = 0; li[0] < n[0]; ++li[0] ) 
                fprintf(stream, "%+11.3e ", ra_get(arr, li));
            fprintf(stream, "\n");    
            break;
        #if RARR_MAX_NDIM > 1
        case 2:
            for ( li[1] = 0; li[1] < n[1]; ++li[1] )
            {
                for ( li[0] = 0; li[0] < n[0]; ++li[0] ) 
                    fprintf(stream, "%+11.3e ", ra_get(arr, li));
                fprintf(stream, "\n");    
            }
            break;
        #endif
        #if RARR_MAX_NDIM > 2
        case 3:
            for ( li[2] = 0; li[2] < n[2]; ++li[2] )
            {
                fprintf(stream, "slice %d\n", li[2]);    
                for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                {
                    for ( li[0] = 0; li[0] < n[0]; ++li[0] ) 
                        fprintf(stream, "%+11.3e ", ra_get(arr, li));
                    fprintf(stream, "\n");    
                }
            }
            break;
        #endif
        default:
            fprintf(stream, "output of arrays of more than 3 dimensions not supported.\n");
            return E_OTHER;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_fwrite(RARR *arr, const char *path)
/*<*
 * Output the working (computational virtual) array to a file
 *
 * Format: binary
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "wb");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = ra_write(arr, stream);
    fclose(stream);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ra_write(RARR *arr, FILE* stream)
/*<*
 * Output the working (computational virtual) array to stream
 *
 * Format: binary
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    IPNT li, n;
    
    ra_size(arr, n);
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            fwrite(arr->_s, sizeof(ireal), n[0], stream);
            break;
        #if RARR_MAX_NDIM > 1
        case 2:
            for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                fwrite(arr->_s + li[1] * arr->_dims[0].n0, sizeof(ireal), n[0], stream);
            break;
        #endif
        #if RARR_MAX_NDIM > 2
        case 3:
            for ( li[2] = 0; li[2] < n[2]; ++li[2] )
            {
                for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                    fwrite(arr->_s + (li[1] + li[2] * arr->_dims[1].n0) * arr->_dims[0].n0, 
                           sizeof(ireal), n[0], stream);
            }
            break;
        #endif
        default:
            return E_OTHER;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
/* D.S. 01.17.10 */

int ra_copy(RARR *arr_des, RARR *arr_src)
/*<*
 * Copy the contents of the source arry to the destination array. 
 * Both arrays must have the same structure and be allocated.
>*/
{ 
  int err = 0;

  if ( (arr_des->_s != NULL) && (arr_src->_s != NULL) ) {
    
    IPNT li, n;
    ra_size(arr_des,n);
    
    switch ( arr_des->ndim ) {
    case 0: 
      break;
      
    case 1:
      memcpy(arr_des->_s, arr_src->_s, n[0]*sizeof(ireal));
      break;
#if RARR_MAX_NDIM > 1
    case 2:
      for ( li[1] = 0; li[1] < n[1]; ++li[1] )
	memcpy(arr_des->_s + li[1]* arr_des->_dims[0].n0, arr_src->_s + li[1]* arr_src->_dims[0].n0, n[0]*sizeof(ireal));
      break;
#endif
#if RARR_MAX_NDIM > 2
    case 3:
      for ( li[2] = 0; li[2] < n[2]; ++li[2] )
	{
	  for ( li[1] = 0; li[1] < n[1]; ++li[1] )
	    memcpy(arr_des->_s + (li[1] + li[2] * arr_des->_dims[1].n0) * arr_des->_dims[0].n0,
		   arr_src->_s + (li[1] + li[2] * arr_src->_dims[1].n0) * arr_src->_dims[0].n0,
		   n[0] * sizeof(ireal));
	}
      break;
#endif
    default:
      return E_OTHER;
    }
    return 0;
  }
  else {
    err = -100;
    return err;
  }
}

/*----------------------------------------------------------------------------*/
/* D.S. 12.17.09 */

int ra_fread(RARR *arr, const char *path)
/*<*
 * read the working (computational virtual) array from a binary file
 *
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "rb");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = ra_read(arr, stream);
    fclose(stream);
    
    return err;
}

/*----------------------------------------------------------------------------*/
/* D.S. 12.12.09 */

int ra_read(RARR *arr, FILE* stream)
/*<*
 * read the working (computational virtual) array from a binary stream
 *
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 >*/
{
    IPNT li, n;
    
    ra_size(arr, n);
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            fread(arr->_s, sizeof(ireal), n[0], stream);
            break;
        #if RARR_MAX_NDIM > 1
        case 2:
            for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                fread(arr->_s + li[1] * arr->_dims[0].n0, sizeof(ireal), n[0], stream);
            break;
        #endif
        #if RARR_MAX_NDIM > 2
        case 3:
            for ( li[2] = 0; li[2] < n[2]; ++li[2] )
            {
                for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                    fread(arr->_s + (li[1] + li[2] * arr->_dims[1].n0) * arr->_dims[0].n0, 
                           sizeof(ireal), n[0], stream);
            }
            break;
        #endif
        default:
            return E_OTHER;
    }
    
    return 0;
}

/*----------------------------------------------------------------------------*/

ireal ra_get(const RARR *arr, IPNT li)
/*<*
 * Get value at a local index relative to gs.
 *
 * [No difference with ra_gget, since gs is alway 0]
 * @param [in] arr - (RARR *)
 * @param [in] li - (IPNT) local index relative to gs
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    int d;
    long pshift = 0L;                   /* pointer shift */
    
    for ( d = arr->ndim - 1; d >= 0; --d ) 
    {
        pshift = pshift * (long)(arr->_dims[d].n0) + (long)(li[d]);

        #ifdef CHECK_BOUNDS
        if ( ra_checkbound(arr, d, li[d]) )
        {
            ra_oob(arr, d, li[d], "ra_get");
            return REAL_NAN;
        }
        #endif
    }
    
    return arr->_s[pshift];
}
/*----------------------------------------------------------------------------*/

ireal ra_gget(const RARR *arr, IPNT gi)
/*<*
 * Get value at a global index.
 *
 * @param [in] arr - (RARR *)
 * @param [in] gi - (IPNT) global index
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    int d;
    const INFODIM *dim;                      /* pointer to current dimension */
    long pshift = 0L;                   /* pointer shift */
    
    for ( d = arr->ndim - 1; d >= 0; --d ) 
    {
        dim = arr->_dims + d;
        pshift = pshift * (long)(dim->n0) + (long)(gi[d] - dim->gs);

        #ifdef CHECK_BOUNDS
        if ( ra_gcheckbound(arr, d, gi[d]) )
        {
            ra_oob(arr, d, gi[d], "ra_gget");
            return REAL_NAN;
        }
        #endif
    }
    
    return arr->_s[pshift];
}
/*----------------------------------------------------------------------------*/

void ra_set(RARR *arr, IPNT li, ireal r)
/*<*
 * Set value at a local index relative to gs.
 *
 * [No difference with ra_gset, since gs is alway 0]
 * @param [out] arr - (RARR *)
 * @param [in] li - (IPNT) local index relative to gs
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    int d;
    long pshift = 0L;                   /* pointer shift */
    
    for ( d = arr->ndim - 1; d >= 0; --d ) 
    {
        pshift = pshift * (long)(arr->_dims[d].n0) + (long)(li[d]);

        #ifdef CHECK_BOUNDS
        if ( ra_checkbound(arr, d, li[d]) )
        {
            ra_oob(arr, d, li[d], "ra_set");
            return;
        }
        #endif
    }
    
    arr->_s[pshift] = r;
}
/*----------------------------------------------------------------------------*/

void ra_gset(RARR *arr, IPNT gi, ireal r)
/*<*
 * Set value at a global index.
 *
 * @param [out] arr - (RARR *)
 * @param [in] gi - (IPNT) global index
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    int d;
    INFODIM *dim;                      /* pointer to current dimension */
    long pshift = 0L;                  /* pointer shift */
    
    for ( d = arr->ndim - 1; d >= 0; --d ) 
    {
        dim = arr->_dims + d;
        pshift = pshift * (long)(dim->n0) + (long)(gi[d] - dim->gs);

        #ifdef CHECK_BOUNDS
        if ( ra_gcheckbound(arr, d, gi[d]) )
        {
            ra_oob(arr, d, gi[d], "ra_gset");
            return;
        }
        #endif
    }
    
    arr->_s[pshift] = r;
}
/*----------------------------------------------------------------------------*/

int ra_size(const RARR *arr, IPNT n)
/*<*
 * Get size of the working (computational virtual) array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] n - (IPNT)
 * @return 0
>*/
{
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) n[d] = arr->_dims[d].n;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_gse(const RARR *arr, IPNT gs, IPNT ge)
/*<*
 * Get the start and end indices of the working (computational virtual) array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] gs - (IPNT) start indices
 * @param [out] ge - (IPNT) end indices
>*/
{
    int d;
    const INFODIM *dim;                      /* pointer to current dimension */
    
    for ( d = 0; d < arr->ndim; ++d )
    {
        dim = arr->_dims + d;
        if ( gs != NULL ) gs[d] = dim->gs;
        if ( ge != NULL ) ge[d] = dim->ge;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_gcheckbound(const RARR *arr, int idim, int gi)
/*<*
 * Check if a gloabal index (and idim) is within bounds of the allocated array.
 *
 * @param [in] arr - (RARR *)
 * @param [in] idim - (int) dimension number. 
 * @param [in] gi - global index
 * @return 0 if within bounds, else error code as in base/include/utils.h.
>*/
{
    const INFODIM *dim;                      /* pointer to current dimension */
    
    if ( (unsigned int)idim >= (unsigned int)(arr->ndim) ) return E_BADDIMINDEX;

    dim = arr->_dims + idim;
    if ( (gi < dim->gs0) || (gi > dim->ge0) ) return E_OUTOFBOUNDS;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_checkbound(const RARR *arr, int idim, int li)
/*<*
 * Check if a local index relative to gs (and idim) is within bounds of the allocated array.
 * 
 * @param [in] arr - (RARR *)
 * @param [in] idim - (int) dimension number.
 * @param [in] li - local index relative to gs
 * @return 0 if within bounds, else error code as in base/include/utils.h.
>*/
{
    int gi;                             /* global index */
    const INFODIM *dim;                      /* pointer to current dimension */
    
    if ( (unsigned int)idim >= (unsigned int)(arr->ndim) ) return E_BADDIMINDEX;

    dim = arr->_dims + idim;
    gi = li + dim->gs;
    if ( (gi < dim->gs0) || (gi > dim->ge0) ) return E_OUTOFBOUNDS;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_ndim(const RARR *arr, int *ndim)
/*<*
 * Get number of dimensions.
 *
 * @param [in] arr - (RARR *)
 * @param [out] ndim - (int *) number of dimensions
 * @return 0
>*/
{
    *ndim = arr->ndim;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_setempty(RARR *arr)
/*< Set the working (computational virtual) array empty >*/
{
    if ( arr->ndim == 0 ) return 0;
	
	return ra_offset_s(arr, IPNT_0, IPNT_0);
}
/*----------------------------------------------------------------------------*/

int ra_empty(RARR *arr, int *empty)
/*<*
 * Array empty query.
 *
 * @param [in] arr - (RARR *)
 * @param [out] empty - (int *) 0: nonempty, 1: empty
 * @return 0
>*/
{
    int d;
    if ( arr->ndim == 0 )
	{
		*empty = 1;
		return 0;
	}
	
    for ( d = 0; d < arr->ndim; ++d ) if ( arr->_dims[d].n <= 0 ) break;
    *empty = (d < arr->ndim);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

#ifdef IWAVE_USE_MPI
/*^*/

int ra_setexchangeinfo(RARR *arr, EXCHANGEINFO *einfo)
/*<*
 * Populates exchange info. 
 *
 * Creates MPI_Datatype inside - do not forget to destoy.
 * @param[in] arr - (RARR *)
 * @param[out] einfo - (EXCHANGEINFO *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
  int ndim;
  IPNT n;
  int err;

  int empty;
	
  einfo->buf = NULL;
  einfo->type = MPI_DATATYPE_NULL;
  ra_empty(arr, &empty);
  if ( empty ) return 0;
  
  
  ra_ndim(arr, &ndim);
  ra_size(arr, n);
  einfo->buf = (void*)arr->_s;
  if ( ndim == 1 )
  {
      err = MPI_Type_contiguous(n[0], IWAVE_MPI_REAL, &(einfo->type));
      if ( err != MPI_SUCCESS ) return E_MPI;
  }
  else if ( ndim == 2 )
  {
      err = MPI_Type_vector(n[1], n[0], arr->_dims[0].n0, IWAVE_MPI_REAL, &(einfo->type));
      if ( err != MPI_SUCCESS ) return E_MPI;
  }
  else if ( ndim == 3 )
  {
      err = MPI_Type_vector(n[1], n[0], arr->_dims[0].n0, IWAVE_MPI_REAL, &(einfo->type2));
      if ( err != MPI_SUCCESS ) return E_MPI;
      err = MPI_Type_hvector(n[2], 1, arr->_dims[0].n0 * arr->_dims[1].n0 * sizeof(ireal), einfo->type2, &(einfo->type));
      if ( err != MPI_SUCCESS )
      {
	  MPI_Type_free(&(einfo->type2));
	  return E_MPI;
      }
  }
  else return E_BADINPUT;
  
  return 0;
}

#endif
/*^*/


/*----------------------------------------------------------------------------*/

int ra_overlap(RARR *arr1, RARR *arr2, int *overlap)
/*<*
 *  Checks if the two working (computational virtual) arrays of arr1 and arr2 overlap.
 * 
 * @param [in] arr1, arr2 - (RARR *)
 * @param [out] overlap - (int *) 0: not overlap, 1: overlap
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
	int e1, e2, idim, ndim1, ndim2;
	IPNT gs1, ge1, gs2, ge2;

	ra_ndim(arr1, &ndim1);
	ra_ndim(arr2, &ndim2);
	if ( ndim1 != ndim2 ) return E_BADINPUT;
	
	*overlap = 0;
	ra_empty(arr1, &e1);
	ra_empty(arr2, &e2);
	if ( e1 || e2 ) return 0; /* empty arrays do not overlap */
	
	ra_gse(arr1, gs1, ge1);
	ra_gse(arr2, gs2, ge2);
	
	for ( idim = 0; idim < ndim1; ++idim )
		if ( (gs1[idim] > ge2[idim]) || (gs2[idim] > ge1[idim]) ) return 0;
		
	*overlap = 1;
	
	return 0;
}
/*----------------------------------------------------------------------------*/

int ra_setoverlap(RARR *arr1, RARR *arr2)
/*<*
 * Set the working (computational virtual) array's dimension info of arr1  to be that of the overlap part of 
 * arr1 and arr2.
 *
 * @param [in,out] arr1 - (RARR *)
 * @param [in] arr2 - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
	int e1, e2, idim, ndim1, ndim2;
	IPNT gs1, ge1, gs2, ge2;

	ra_ndim(arr1, &ndim1);
	ra_ndim(arr2, &ndim2);
	if ( ndim1 != ndim2 ) return E_BADINPUT;
	
	ra_empty(arr1, &e1);
	ra_empty(arr2, &e2);
	if ( e1 || e2 ) return ra_setempty(arr1); /* empty array => empty overlap */
	
	ra_gse(arr1, gs1, ge1);
	ra_gse(arr2, gs2, ge2);
	
	for ( idim = 0; idim < ndim1; ++idim )
	{
		gs1[idim] = iwave_max( gs1[idim], gs2[idim] );
		ge1[idim] = iwave_min( ge1[idim], ge2[idim] );
		if ( gs1[idim] > ge1[idim] ) return ra_setempty(arr1);
	}
		
	return ra_greset(arr1, gs1, ge1);
}
/*----------------------------------------------------------------------------*/

int ra_zero(RARR *arr)
/*< Set the entries of the allocated arry all zero >*/
{
	int len, dim, i;

	ra_ndim(arr, &dim);
	if ( dim <= 0 ) return 0;

	len = 1;
	for ( i=0; i < dim; ++i ) len *= arr->_dims[i].n0;
	for ( i=0; i < len; ++i ) arr->_s0[i] = REAL_ZERO;

	return 0;
}
/*----------------------------------------------------------------------------*/

int ra_deepcopy(RARR const * src, RARR * tgt)
/*<*
 * Copy all of the members of the source struct to the target struct 
 *
 * No allocation occurs. tgt works as a reference of src.
 * 
 * @param [in] src - (RARR const *) src cannot be pointed to other array, but one 
 *                    can change the array pointed by it
 * @param [out] tgt - (RARR *)
 * @return 0
>*/
{
	int i;
        tgt->ndim = src->ndim;
        tgt->_s = src->_s;
        tgt->_s0 = src->_s0; 
        for (i=0; i<src->ndim; ++i)
	{       
	    tgt->_dims[i].n = src->_dims[i].n;
	    tgt->_dims[i].n0 = src->_dims[i].n0;
            tgt->_dims[i].gs = src->_dims[i].gs;
	    tgt->_dims[i].gs0 = src->_dims[i].gs0;
	    tgt->_dims[i].ge = src->_dims[i].ge;
            tgt->_dims[i].ge0 = src->_dims[i].ge0;
	}
    
	return 0;
}
