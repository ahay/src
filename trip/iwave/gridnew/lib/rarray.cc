/* 
rarray.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "rarray.h"

/* #undef VERBOSE */
#define VERBOSE

/*----------------------------------------------------------------------------*/
/*
Out of bounds message.
*/
static void ra_oob(const RARR *arr, int idim, int ind, const char *funname)
{
#ifdef VERBOSE
    fprintf(stderr, "%s OUT OF BOUNDS: dim = %d, ind = %d\n", funname, idim, ind);
    ra_dump(arr, stderr);
#endif
}
/*----------------------------------------------------------------------------*/

int ra_getdim(IPNT s, IPNT e) {
  int ndim = 0;
  for (int d = RARR_MAX_NDIM-1; d > -1; d--) 
    if ((ndim==0)  && (e[d]-s[d] > 0)) ndim=d+1;
  return ndim;
}

/*----------------------------------------------------------------------------*/

int ra_setnull(RARR *arr)
{
    memset((void*)arr, 0, sizeof(RARR));
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_create(RARR *arr, IPNT gs, IPNT ge)
{
    int err;
    
    err = ra_declare(arr, gs, ge);
    if ( err ) {    
      fprintf(stderr,"ERROR: RARR::ra_create - from ra_declare err=%d\n",err);
      ra_setnull(arr);
      return err;
    }
    err = ra_allocate(arr);

    if ( err )
    {    
      fprintf(stderr,"ERROR: RARR::ra_create - from ra_allocate err=%d\n",err);
        ra_setnull(arr);
        return err;
    }
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_create_s(RARR *arr, IPNT gs, IPNT n)
{
    IPNT ge; 
    int d;
    
    for ( d = 0; d < RARR_MAX_NDIM; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_create(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_create_e(RARR *arr, IPNT ge, IPNT n)
{
    IPNT gs;
    int d;
    
    for ( d = 0; d < RARR_MAX_NDIM; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_create(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare_s(RARR *arr, IPNT gs, IPNT n)
{
    IPNT ge; 
    int d;
    
    for ( d = 0; d < RARR_MAX_NDIM; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_declare(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare_e(RARR *arr, IPNT ge, IPNT n)
{
    IPNT gs;
    int d;
    
    for ( d = 0; d < RARR_MAX_NDIM; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_declare(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_declare(RARR *arr, IPNT gs, IPNT ge)
{
    int d;
    INFODIM *dim;                      /* pointer to current dimension */
    
    ra_setnull(arr);                   /* empty array */

    for ( d = 0; d < RARR_MAX_NDIM; ++d )
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
    
    arr->ndim = ra_getdim(gs,ge);

    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_allocate(RARR *arr) {
  int d;                              /* counter */
  long size;                          /* size to allocate */
#if RARR_MAX_NDIM > 1
  long psize;                         /* size of slice    */
  IPNT s0;                            /* start indices    */
  IPNT e0;                            /* end indices      */
  IPNT n0;                            /* axis lengths     */
#endif
  
  if ( arr->_s0 != NULL ) return E_ALREADYALLOC;
  
  size = 1L;                          /* allocation size */
  for ( d = 0; d < arr->ndim; ++d ) size *= (long)(arr->_dims[d].n0);
  
  if ( size > 0L ) {                   /* allocate memory */

    ra_a_gse(arr,s0,e0);
    ra_a_size(arr,n0);

#if RARR_MAX_NDIM > 0 

    arr->_s0 = (ireal*)usermalloc_(size * sizeof(ireal));
    if ( arr->_s0 == NULL ) return E_ALLOC;
    arr->_s = arr->_s0;
    arr->_s1 = arr->_s0-s0[0];
    
#endif
    
    /* added 03.11.12 WWS: multidim array access */
    
#if RARR_MAX_NDIM > 1

    psize = 1L;

    for (d=1; d<arr->ndim; ++d) psize *= (long)(arr->_dims[d].n0);
    if (psize > 0L && arr->ndim > 1) {
      arr->_s02 = (ireal **)usermalloc_(psize * sizeof(ireal*));
      arr->_s2 = arr->_s02-s0[1];
      for (d=0;d<psize; ++d)  {
	arr->_s02[d]=&(arr->_s0[d*n0[0]]);
	arr->_s2[d+s0[1]]=arr->_s02[d]-s0[0];
      }
    }
    else {
      arr->_s02=NULL;
      arr->_s2=NULL;
    }
#endif

#if RARR_MAX_NDIM > 2

    psize = 1L;
    for (d=2; d<arr->ndim; ++d) psize *= (long)(arr->_dims[d].n0);
    if ((psize > 0L) && arr->ndim > 2) {
      arr->_s03 = (ireal ***)usermalloc_(psize * sizeof(ireal**));
      arr->_s3 = arr->_s03-s0[2];
      for (d=0;d<psize; ++d) {
	arr->_s03[d]=&(arr->_s02[d*n0[1]]);
	arr->_s3[d+s0[2]]=arr->_s03[d]-s0[1];
      }
    }
    else {
      arr->_s03=NULL;
      arr->_s3=NULL;
    }
#endif

#if RARR_MAX_NDIM > 3

    psize = 1L;
    for (d=3; d<arr->ndim; ++d) psize *= (long)(arr->_dims[d].n0);
    if (psize > 0L && arr->ndim > 3) {
      arr->_s04 = (ireal ****)usermalloc_(psize * sizeof(ireal***));
      arr->_s4 = arr->_s04-s0[3];
      for (d=0;d<psize; ++d) {
	arr->_s04[d]=&(arr->_s03[d*n0[2]]);
	arr->_s4[d+s0[3]]=arr->_s04[d]-s0[2];
      }
    }
    else {
      arr->_s04=NULL;
      arr->_s4=NULL;
    }

#endif

#if RARR_MAX_NDIM > 4

    psize = 1L;
    for (d=4; d<arr->ndim; ++d) psize *= (long)(arr->_dims[d].n0);
    if (psize > 0L && arr->ndim > 4) {
      arr->_s05 = (ireal *****)usermalloc_(psize * sizeof(ireal****));
      arr->_s5 = arr->_s05-s0[4];
      for (d=0;d<psize; ++d) {
	arr->_s05[d]=&(arr->_s04[d*n0[3]]);
	arr->_s5[d+s0[4]]=arr->_s05[d]-s0[3];
      }
    }
    else {
      arr->_s04=NULL;
      arr->_s4=NULL;
    }

#endif

  }

  return 0;
}
/*----------------------------------------------------------------------------*/

int ra_destroy(RARR *arr) {
  if (arr->_s0) userfree_(arr->_s0);
#if RARR_MAX_NDIM > 1
  if (arr->_s02) userfree_(arr->_s02);
#endif
#if RARR_MAX_NDIM > 2
  if (arr->_s03) userfree_(arr->_s03);
#endif
#if RARR_MAX_NDIM > 3
  if (arr->_s04) userfree_(arr->_s04);
#endif
#if RARR_MAX_NDIM > 4
  if (arr->_s05) userfree_(arr->_s05);
#endif

  ra_setnull(arr);

  return 0;
}
/*----------------------------------------------------------------------------*/

int ra_offset(RARR *arr, const IPNT os, const IPNT oe)
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

int ra_greset(RARR *arr, const IPNT gs, const IPNT ge) {

  int d;
  /*  IPNT idx;                         index workspace */
  INFODIM *dim;                      /* pointer to current dimension */
  long pshift=0L;                    /* pointer shift */
  long soff[RARR_MAX_NDIM+1];        /* pointer offset - start */
  /* long eoff[RARR_MAX_NDIM+1];        pointer offset - end */
  RARR arrold = *arr;                /* remember old to restore if error */
    
  /* cycle backwards to compute pshift correctly */
  for ( d = arr->ndim - 1; d >= 0; --d ) {
    dim = arr->_dims + d;
    dim->gs = gs[d];
    dim->ge = ge[d];
    dim->n = dim->ge - dim->gs + 1;
    if ( (dim->n < 0) || (dim->gs < dim->gs0) || (dim->ge > dim->ge0) ) {
      fprintf(stderr,"Error: ra_greset - out of bounds\n");
      fprintf(stderr,"n=%d gs=%d gs0=%d ge=%d ge0=%d\n",dim->n,dim->gs,dim->gs0,dim->ge,dim->ge0);
      *arr = arrold;
      return E_OUTOFBOUNDS;
    }
    soff[d] = (long)(dim->gs - dim->gs0);
    /*
    eoff[d] = (long)(dim->ge0 - dim->ge);
    */
    pshift = pshift * (long)(dim->n0) + soff[d];

  }
  
  if ( arr->_s != NULL ) arr->_s = arr->_s0 + pshift;

  return 0;
}
/*----------------------------------------------------------------------------*/

int ra_greset_s(RARR *arr, const IPNT gs, const IPNT n)
{
    IPNT ge;
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) ge[d] = gs[d] + n[d] - 1;
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_greset_e(RARR *arr, const IPNT ge, const IPNT n)
{
    IPNT gs; 
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) gs[d] = ge[d] - n[d] + 1;
    
    return ra_greset(arr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int ra_dump(const RARR *arr, FILE* stream)
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
/* WWS. 14.02.12 */
int ra_a_copy(RARR *arr_des, RARR *arr_src)
{ 
  int err = 0;
  int i, n;
  IPNT n_des, n_src;

  if ( !(arr_des->_s) ) {
    fprintf(stderr,"ERROR: RARR::ra_a_copy\n");
    fprintf(stderr,"  destination array not initialized\n");
  }
  if ( !(arr_des->_s) ) {
    fprintf(stderr,"ERROR: RARR::ra_a_copy\n");
    fprintf(stderr,"  destination array not initialized\n");
  }
  if ( !(arr_src->_s) ) {
    fprintf(stderr,"ERROR: RARR::ra_a_copy\n");
    fprintf(stderr,"  source array not initialized\n");
  }
  if ( (arr_des->ndim != arr_src->ndim) ) {
    fprintf(stderr,"ERROR: RARR::ra_a_copy\n");
    fprintf(stderr,"  dest ndim = %d != src ndim = %d\n",arr_des->ndim, arr_src->ndim);
  }
  if ( (arr_des->_s != NULL) && (arr_src->_s != NULL) &&
       (arr_des->ndim == arr_src->ndim) ) {
    
    ra_a_size(arr_des,n_des);
    ra_a_size(arr_src,n_src);

    n=1;
    for (i=0;i<arr_des->ndim;i++) {
      if (n_des[i]==n_src[i]) n *= n_des[i];
      else {
	fprintf(stderr,"ERROR: RARR::ra_a_copy\n");
	fprintf(stderr,"  dest n[%d]=%d != src n[%d]=%d\n",i,n_des[i],i,n_src[i]);
	return E_BADINDEX;
      }
    }
    memcpy(arr_des->_s0, arr_src->_s0, n*sizeof(ireal));
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
{
    IPNT li, n;
    
    ra_size(arr, n);
    
    switch ( arr->ndim ) 
    {
        case 0: 
            break;
        
        case 1:
            if ((size_t) n[0] != fread(arr->_s, sizeof(ireal), n[0], stream))
		return E_OTHER;
            break;
        #if RARR_MAX_NDIM > 1
        case 2:
            for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                if ((size_t) n[0] != fread(arr->_s + li[1] * arr->_dims[0].n0, sizeof(ireal), n[0], stream))
		    return E_OTHER;
            break;
        #endif
        #if RARR_MAX_NDIM > 2
        case 3:
            for ( li[2] = 0; li[2] < n[2]; ++li[2] )
            {
                for ( li[1] = 0; li[1] < n[1]; ++li[1] )
                    if ((size_t) n[0] != fread(arr->_s + (li[1] + li[2] * arr->_dims[1].n0) * arr->_dims[0].n0, sizeof(ireal), n[0], stream))
			return E_OTHER;
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
{
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) n[d] = arr->_dims[d].n;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_datasize(const RARR *arr, size_t * n)
{
    int d;
    *n = 1;
    for ( d = 0; d < arr->ndim; ++d ) *n *= arr->_dims[d].n;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_a_size(const RARR *arr, IPNT n)
{
    int d;
    
    for ( d = 0; d < arr->ndim; ++d ) n[d] = arr->_dims[d].n0;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_a_datasize(const RARR *arr, size_t * n)
{
    int d;
    *n=1;
    for ( d = 0; d < arr->ndim; ++d ) *n *= arr->_dims[d].n0;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_gse(const RARR *arr, IPNT gs, IPNT ge)
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

int ra_se(const RARR *arr, IPNT s, IPNT e)
{
    int d;
    const INFODIM *dim;                      /* pointer to current dimension */
    
    for ( d = 0; d < arr->ndim; ++d )
    {
        dim = arr->_dims + d;
        if ( s != NULL ) s[d] = dim->gs-dim->gs0;
        if ( e != NULL ) e[d] = dim->ge-dim->gs0;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_a_gse(const RARR *arr, IPNT gs, IPNT ge)
{
    int d;
    const INFODIM *dim;                      /* pointer to current dimension */
    
    IASN(gs,IPNT_1);
    IASN(ge,IPNT_0);

    for ( d = 0; d < arr->ndim; ++d )
    {
        dim = arr->_dims + d;
        if ( gs != NULL ) gs[d] = dim->gs0;
        if ( ge != NULL ) ge[d] = dim->ge0;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_ds(const RARR * tgt, const RARR * src, IPNT ds) {

  int d;
  const INFODIM *tdim;                      /* pointer to target dimension */
  const INFODIM *sdim;                      /* pointer to source dimension */

  if (tgt->ndim != src->ndim) return E_BADINPUT;

  for (d=0; d<tgt->ndim; d++) {
    tdim = tgt->_dims + d;
    sdim = src->_dims + d;
    ds[d] = sdim->gs0 - tdim->gs0;      
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
int ra_gcheckbound(const RARR *arr, int idim, int gi)
{
    const INFODIM *dim;                      /* pointer to current dimension */
    
    if ( (unsigned int)idim >= (unsigned int)(arr->ndim) ) return E_BADDIMINDEX;

    dim = arr->_dims + idim;
    if ( (gi < dim->gs0) || (gi > dim->ge0) ) return E_OUTOFBOUNDS;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_checkbound(const RARR *arr, int idim, int li)
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
{
    *ndim = arr->ndim;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ra_setempty(RARR *arr)
{
    if ( arr->ndim == 0 ) return 0;
	
	return ra_offset_s(arr, IPNT_0, IPNT_0);
}
/*----------------------------------------------------------------------------*/

int ra_empty(RARR *arr, int *empty)
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

int ra_setexchangeinfo(RARR *arr, EXCHANGEINFO *einfo)
{
#ifdef IWAVE_USE_MPI
  int ndim;
  IPNT n;
  int err;
#endif

  int empty;
	
	einfo->buf = NULL;
	einfo->type = MPI_DATATYPE_NULL;
	ra_empty(arr, &empty);
	if ( empty ) return 0;

	#ifdef IWAVE_USE_MPI
	
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
	
	#endif
	
	return 0;
}
/*----------------------------------------------------------------------------*/

int ra_overlap(RARR *arr1, RARR *arr2, int *overlap)
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

int ra_zero(RARR * a) {

    IPNT li;         /* counter */
    IPNT n,n0;   /* sizes */
  
  if (a->_s != NULL) {
    
    ra_size(a,n);
    ra_a_size(a,n0);

    switch ( a->ndim ) {
    case 0: 
      break;
      
    case 1:
      for (li[0]=0;li[0]<n[0];li[0]++) (a->_s)[li[0]] = REAL_ZERO;
      break;
#if RARR_MAX_NDIM > 1
    case 2:
      for ( li[1] = 0; li[1] < n[1]; ++li[1] )
	for (li[0]=0;li[0]<n[0];li[0]++) (a->_s)[li[0]+li[1]*n0[0]] = REAL_ZERO;
      break;
#endif
#if RARR_MAX_NDIM > 2
    case 3:
      for ( li[2] = 0; li[2] < n[2]; ++li[2] ) {
	for ( li[1] = 0; li[1] < n[1]; ++li[1] )
	  for (li[0]=0;li[0]<n[0];li[0]++) (a->_s)[li[0]+(li[1]+li[2]*n0[1])*n0[0]] = REAL_ZERO;
      }
      break;
#endif
    default:
      return E_OTHER;
    }
    
    return 0;
  }
  return E_BADINDEX;
}
      
/*----------------------------------------------------------------------------*/

int ra_a_zero(RARR *arr)
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

int ra_a_scale(RARR *arr, ireal fac)
{
	int len, dim, i;

	ra_ndim(arr, &dim);
	if ( dim <= 0 ) return 0;

	len = 1;
	for ( i=0; i < dim; ++i ) len *= arr->_dims[i].n0;
	for ( i=0; i < len; ++i ) arr->_s0[i] *= fac;

	return 0;
}
/*----------------------------------------------------------------------------*/
int ra_deepcopy(RARR const * src, RARR * tgt)
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

int ra_a_inner(RARR const * arr1, RARR const * arr2, ireal * ip) {

    int i;         /* counter */
    int n;         /* total data length */
    IPNT n1, n2;   /* sizes */
  
  if ( (arr2->_s0 != NULL) && (arr2->_s0 != NULL) &&
       (arr1->ndim == arr2->ndim) ) {
    
    ra_a_size(arr1,n1);
    ra_a_size(arr2,n2);
    n=1;
    for (i=0;i<arr1->ndim;i++) {
      if (n1[i]==n2[i]) n *= n1[i];
      else return E_BADINDEX;
    }

    *ip=REAL_ZERO;
    for (i=0;i<n;i++) *ip+=(arr2->_s0)[i]*(arr1->_s0)[i];
    return 0;
  }
  return E_BADINDEX;
}
      
int ra_axpy(RARR * arry, RARR const * arrx, ireal a) {

  int i;
  IPNT li;         /* counter */
  IPNT nx, ny;   /* sizes */
  IPNT nx0, ny0; /* allocated sizes */
  
  if ( (arrx->_s != NULL) && (arry->_s != NULL) &&
       (arrx->ndim == arry->ndim) ) {
    
    ra_size(arrx,nx);
    ra_size(arry,ny);
    ra_a_size(arrx,nx0);
    ra_a_size(arry,ny0);

    for (i=0;i<arrx->ndim;i++) 
      if (nx[i]!=ny[i]) return E_BADINDEX;

    switch ( arrx->ndim ) {
    case 0: 
      break;
      
    case 1:
      for (li[0]=0;li[0]<nx[0];li[0]++) (arry->_s)[li[0]] += a*((arrx->_s)[li[0]]);            
      break;
#if RARR_MAX_NDIM > 1
    case 2:
      for ( li[1] = 0; li[1] < nx[1]; ++li[1] )
	for (li[0]=0;li[0]<nx[0];li[0]++) (arry->_s)[li[0]+li[1]*ny0[0]] += a*((arrx->_s)[li[0]+li[1]*nx0[0]]);            
      break;
#endif
#if RARR_MAX_NDIM > 2
    case 3:
      for ( li[2] = 0; li[2] < nx[2]; ++li[2] ) {
	for ( li[1] = 0; li[1] < nx[1]; ++li[1] )
	  for (li[0]=0;li[0]<nx[0];li[0]++) (arry->_s)[li[0]+(li[1]+li[2]*ny0[1])*ny0[0]] += a*((arrx->_s)[li[0]+(li[1]+li[2]*nx0[1])*nx0[0]]);            
      }
      break;
#endif
    default:
      return E_OTHER;
    }
    
    return 0;
  }
  return E_BADINDEX;
}
      
int ra_compare_meta(const RARR * a, const RARR * b) {

  int i;
  int ndima;
  int ndimb;
  IPNT gea;
  IPNT gsa;
  IPNT geb;
  IPNT gsb;
  
  ra_ndim(a,&ndima);
  ra_ndim(b,&ndimb);

  if (ndima != ndimb) return 1;

  ra_a_gse(a,gsa,gea);
  ra_a_gse(b,gsb,geb);

  for (i=0;i<ndima;i++) {
    if (gsa[i] != gsb[i]) return 2+2*i;
    if (gea[i] != geb[i]) return 2+2*i+1;
  }
  return 0
;

}
