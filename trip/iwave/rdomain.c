/* Domain type and its operations. */
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
rdomain.c
Igor Terentyev.
*/
/*============================================================================*/

#include <trip/base.h>
/*^*/

#include "rdomain.h"

#include "rarray.h"
#include "exchangeinfo.h"
/*^*/

#ifndef _sf_rdomain_h

typedef struct
{
  /** number of arrays in the domain */
  int narr;
  /** array storage pointer. RDOM_MAX_NARR: 
   * maximum number of arrays in one domain.*/
  RARR _s[RDOM_MAX_NARR];
} RDOM;
/*^*/

#endif

/*----------------------------------------------------------------------------*/
/*
Bad array index function.
*/
static void rd_bai(const RDOM *dom, int iarr, const char *funname)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) )
    {
        fprintf(stderr, "%s BAD ARRAY INDEX: ind = %d\n", funname, iarr);
        rd_a_dump(dom, stderr);
    }
}
/*----------------------------------------------------------------------------*/
/*
Creates/declares domain using fun to create/declare arrays.
*/
static int rd_a_create_(RDOM *dom, int narr, int ndim, IPNT v1[], IPNT v2[], RA_CREATE_FUN fun)
{
    int a, err;
    //    fprintf(stderr,"rd_a_create_: narr=%d\n",narr);

    rd_a_setnull(dom);
    
    if ( (unsigned int)narr > (unsigned int)RDOM_MAX_NARR ) return E_BADINPUT;
    dom->narr = narr;
    
    for ( a = 0; a < narr; ++a )
    {
        err = fun(dom->_s + a, ndim, v1[a], v2[a]);
        if ( err )
        {
            rd_a_destroy(dom);
            return err;
        }
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
/*
Creates/declares and adds array to domain using fun to create/declare array.
*/
static int rd_create_(RDOM *dom, int ndim, IPNT v1, IPNT v2, RA_CREATE_FUN fun)
{
    int err;
    
    if ( dom->narr >= RDOM_MAX_NARR ) return E_BADARRINDEX;
    
    err = fun(dom->_s + dom->narr, ndim, v1, v2);
    if ( err ) return err;

    ++(dom->narr);
    
    return 0;
}
/*----------------------------------------------------------------------------*/
/*
Checks array index and calls array set function.
*/
static int rd_set_(RDOM *dom, int iarr, IPNT v1, IPNT v2, RA_SET_FUN fun)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return fun(dom->_s + iarr, v1, v2);
}
/*----------------------------------------------------------------------------*/

int rd_a_setnull(RDOM *dom)
/*< Set domain (all fields) to zeros. >*/
{
    memset((void*)dom, 0, sizeof(RDOM));    

    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_setnull(RDOM *dom, int iarr)
/*<*
 * Set all fields in the specified array of the domain to zeros
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in]  iarr - (int) array index
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h. One of the errors arises when iarr>=narr.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setnull(dom->_s + iarr);
}
/*----------------------------------------------------------------------------*/

int rd_a_create_s(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])
/*<*
 * Create all arrays in a domain (given start indices (dgs) and sizes (dn)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_a_create_(dom, narr, ndim, dgs, dn, ra_create_s);
}
/*----------------------------------------------------------------------------*/

int rd_a_create_e(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])
/*<*
 * Create all arrays in a domain (given end indices (dge) and sizes (dn)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_a_create_(dom, narr, ndim, dge, dn, ra_create_e);
}
/*----------------------------------------------------------------------------*/

int rd_a_create(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])
/*<*
 * Create all arrays in a domain (given start indices (dgs) and end indices (dge)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
  //  fprintf(stderr,"rd_a_create\n");
    return rd_a_create_(dom, narr, ndim, dgs, dge, ra_create);
}
/*----------------------------------------------------------------------------*/

int rd_create_s(RDOM *dom, int ndim, IPNT gs, IPNT n)
/*<*
 * Create next array in a domain (given start indices (gs) and sizes (n) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_create_(dom, ndim, gs, n, ra_create_s);
}
/*----------------------------------------------------------------------------*/

int rd_create_e(RDOM *dom, int ndim, IPNT ge, IPNT n)
/*<*
 * Create next array in a domain (given end indices (ge) and sizes (n) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 >*/
{
    return rd_create_(dom, ndim, ge, n, ra_create_e);
}
/*----------------------------------------------------------------------------*/

int rd_create(RDOM *dom, int ndim, IPNT gs, IPNT ge)
/*<*
 * Create next array in a domain (given start indices (gs) and end indices (ge) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_create_(dom, ndim, gs, ge, ra_create);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare_s(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])
/*<*
 * Declare all arrays in a domain (given dgs and dn).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_a_create_(dom, narr, ndim, dgs, dn, ra_declare_s);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare_e(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])
/*<*
 * Declare all arrays in a domain (given dge and dn).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_a_create_(dom, narr, ndim, dge, dn, ra_declare_e);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])
/*<*
 * Declare all arrays in a domain (given dgs and dge).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] ndim - (int) number of dimensions
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_a_create_(dom, narr, ndim, dgs, dge, ra_declare);
}
/*----------------------------------------------------------------------------*/

int rd_declare_s(RDOM *dom, int ndim, IPNT gs, IPNT n)
/*<*
 * Declare the next array in a domain (given gs and n).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_create_(dom, ndim, gs, n, ra_declare_s);
}
/*----------------------------------------------------------------------------*/

int rd_declare_e(RDOM *dom, int ndim, IPNT ge, IPNT n)
/*<*
 * Declare the next array in a domain (given ge and n).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_create_(dom, ndim, ge, n, ra_declare_e);
}
/*----------------------------------------------------------------------------*/

int rd_declare(RDOM *dom, int ndim, IPNT gs, IPNT ge)
/*<*
 * Declare the next array in a domain (given gs and ge).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_create_(dom, ndim, gs, ge, ra_declare);
}

/*----------------------------------------------------------------------------*/
/*
Allocates all arrays.
*/

int rd_a_allocate(RDOM *dom)
/*< Allocate memory for all arrays in a domain. >*/
{
    int a, err;
    
    for ( a = 0; a < dom->narr; ++a )
    {
        err = ra_allocate(dom->_s + a);
        if ( err && (err != E_ALREADYALLOC) ) /* supress E_ALREADYALLOC error */
        {
            rd_a_destroy(dom);
            return err;
        }
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
/*
Allocates array.
*/

int rd_allocate(RDOM *dom, int iarr)
/*<* 
 * Allocate memory for a specified array in a domain.
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h. One of the errors arises when iarr>=narr.
>*/
{
    int err;
    
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    err = ra_allocate(dom->_s + iarr);
    if ( err ) return err;

    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_a_destroy(RDOM *dom)
/*< Destroy domain (STORAGE DEALLOCATION). >*/
{
    int a;
    
    for ( a = 0; a < dom->narr; ++a ) {
      /*
      fprintf(stderr,"rd_a_destroy: a=%d\n",a);
      */
      ra_destroy(dom->_s + a);
    }
    rd_a_setnull(dom);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_greset_s(RDOM *dom, int iarr, IPNT gs, IPNT n)
/*< Reset a specified working (computational virtual) arrays in a domain (given gs and n) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_set_(dom, iarr, gs, n, ra_greset_s);
}
/*----------------------------------------------------------------------------*/

int rd_greset_e(RDOM *dom, int iarr, IPNT ge, IPNT n)
/*<*
 * Reset a specified working (computational virtual) arrays in a domain (given ge and n) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    return rd_set_(dom, iarr, ge, n, ra_greset_e);
}
/*----------------------------------------------------------------------------*/

int rd_greset(RDOM *dom, int iarr, IPNT gs, IPNT ge)
/*<*
 * Reset a specified working (computational virtual) array in a domain (given os and n) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset_s.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] os  - (IPNT) start index offsets (forward) of the working (computational virtual) array
 * @param [in] 	n  - (IPNT) sizes of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
>*/
{
    return rd_set_(dom, iarr, gs, ge, ra_greset);
}
/*----------------------------------------------------------------------------*/

int rd_offset_s(RDOM *dom, int iarr, IPNT os, IPNT n)
/*<*
 * Reset a specified working (computational virtual) array in a domain (given os and n) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset_s.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] os  - (IPNT) start index offsets (forward) of the working (computational virtual) array
 * @param [in] 	n  - (IPNT) sizes of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 >*/
{
    return rd_set_(dom, iarr, os, n, ra_offset_s);
}
/*----------------------------------------------------------------------------*/

int rd_offset_e(RDOM *dom, int iarr, IPNT oe, IPNT n)
/*<*
 * Reset a specified working (computational virtual) array in a domain (given oe and n) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset_e.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] oe - (IPNT) end index offsets (backward) of the working (computational virtual) array
 * @param [in] 	n  - (IPNT) sizes of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
>*/
{
    return rd_set_(dom, iarr, oe, n, ra_offset_e);
}
/*----------------------------------------------------------------------------*/

int rd_offset(RDOM *dom, int iarr, IPNT os, IPNT oe)
/*<*
 * Reset a specified working (computational virtual) array in a domain (given os and oe) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] os  - (IPNT) start index offsets (forward) of the working (computational virtual) array
 * @param [in] oe - (IPNT) end index offsets (backward) of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
>*/
{
    return rd_set_(dom, iarr, os, oe, ra_offset);
}
/*----------------------------------------------------------------------------*/

int rd_a_greset(RDOM *dom, IPNT dgs[], IPNT dge[])
/*<*
 * Reset all the working (computational virtual) arrays in a domain (given dgs and dge) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every working (computational virtual) array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every working (computational virtual) array
 * * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
>*/
{
    int a, err;
    
    for ( a = 0; a < dom->narr; ++a )
    {
        err = ra_greset(dom->_s + a, dgs[a], dge[a]);
        if ( err ) return err;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_a_dump(const RDOM *dom, FILE* stream)
/*< Dump information of all arrays in a domain  >*/
{
    int a;
    
    fprintf(stream, "domain dump: narr = %d\n", dom->narr);

    for ( a = 0; a < dom->narr; ++a )
    {
        fprintf(stream, "arr %d\n", a);
        ra_dump(dom->_s + a, stream);
        fprintf(stream, "----------------------------------------\n");
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_dump(const RDOM *dom, int iarr, FILE* stream)
/*<*
 * Dump information of a specified arrays in a domain
 *
 * @param [in] dom - (const RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_dump(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_a_print(RDOM *dom, FILE* stream)
/*<*
 * Output all the working (computational virtual) arrays in a domain to a stream
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    int a, err;
    
    for ( a = 0; a < dom->narr; ++a )
    {
        fprintf(stream, "arr %d\n", a);
        err = ra_print(dom->_s + a, stream);
        if ( err ) return err;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_a_fprint(RDOM *dom, const char *path)
/*<*
 * Output all the working (computational virtual) arrays in a domain to a file
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    FILE *stream;
    int err;

    stream = fopen(path, "w");
    if ( stream == NULL ) return E_FILEOPEN;
        
    err = rd_a_print(dom, stream);
    fclose(stream);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int rd_print(RDOM *dom, int iarr, FILE* stream)
/*<*
 * Output a specified working (computational virtual) array in a domain to a stream
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_print(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_fprint(RDOM *dom, int iarr, const char *path)
/*<*
 * Output a specified working (computational virtual) arrays in a domain to a file
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fprint(dom->_s + iarr, path);
}
/*----------------------------------------------------------------------------*/

int rd_write(RDOM *dom, int iarr, FILE* stream)
/*<*
 * Output a specified working (computational virtual) array in a domain to a stream
 *
 * Format: binary
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_write(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_fwrite(RDOM *dom, int iarr, const char *path)
/*<*
 * Output a specified working (computational virtual) array in a domain to a file
 *
 * Format: binary
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fwrite(dom->_s + iarr, path);
}
/*----------------------------------------------------------------------------*/

int rd_printslice(RDOM *dom, int iarr, FILE* stream, int idim, int li)
/*<*
 * Output a slice of a specified working (computational virtual) array in a domain to a stream.
 *
 * check if iarr < narr. Then call \ref ra_printslice.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_printslice(dom->_s + iarr, stream, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_fprintslice(RDOM *dom, int iarr, const char *path, int idim, int li)
/*<*
 * Output a slice of a specified working (computational virtual) array in a domain to a stream.
 *
 * check if iarr < narr. Then call \ref ra_fprintslice.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fprintslice(dom->_s + iarr, path, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_writeslice(RDOM *dom, int iarr, FILE* stream, int idim, int li)
/*<*
 * Output a slice of a specified working (computational virtual) array in a domain to a binary stream.
 *
 * check if iarr < narr. Then call \ref ra_writeslice.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_writeslice(dom->_s + iarr, stream, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_fwriteslice(RDOM *dom, int iarr, const char *path, int idim, int li)
/*<*
 * Output a slice of a specified working (computational virtual) array in a domain to a binary file.
 *
 * check if iarr < narr. Then call \ref ra_fwriteslice.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fwriteslice(dom->_s + iarr, path, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_a_fsprint(RDOM *dom, const char *path)
/*<*
 * Output each the working (computational virtual) array in a domain to a corresponding file
 *
 * the i'th array is stored in the file named 'str(path)+str(i)'
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] path - (const char *) file name (general)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 >*/
{
    int a, err;
    char path2[70];
    
    if ( strlen(path) > 60L ) return E_BADINPUT;
        
    for ( a = 0; a < dom->narr; ++a )
    {
        if ( sprintf(path2, "%s%d", path, a) < 0 ) return E_OTHER;
            
        err = ra_fprint(dom->_s + a, path2);
        if ( err ) return err;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

ireal rd_get(const RDOM *dom, int iarr, IPNT li)
/*<*
 * Get value at a local index relative to gs in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_get.
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] li - (IPNT) local index relative to gs
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_get");
    #endif
    
    return ra_get(dom->_s + iarr, li);
}
/*----------------------------------------------------------------------------*/

ireal rd_gget(const RDOM *dom, int iarr, IPNT gi)
/*<*
 * Get value at a global index in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_gget.
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gi - (IPNT) global index
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_gget");
    #endif
    
    return ra_gget(dom->_s + iarr, gi);
}
/*----------------------------------------------------------------------------*/

void rd_set(RDOM *dom, int iarr, IPNT li, ireal r)
/*<*
 * Set value at a local index relative to gs in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_set.
 * [No difference with ra_gset, since gs is alway 0]
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] li - (IPNT) local index relative to gs
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_set");
    #endif
    
    ra_set(dom->_s + iarr, li, r);
}
/*----------------------------------------------------------------------------*/

void rd_gset(RDOM *dom, int iarr, IPNT gi, ireal r)
/*<*
 * Set value at a global index in a specified working (computational virtual) array
 *
 * Refer to \ref ra_gset.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gi - (IPNT) global index
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
>*/
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_gset");
    #endif
    
    ra_gset(dom->_s + iarr, gi, r);
}
/*----------------------------------------------------------------------------*/

int rd_size(RDOM *dom, int iarr, IPNT n)
/*<*
 * Get size of a specified working (computational virtual) arry in a domain.
 *
 * Refer to \ref ra_size. 
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_size(dom->_s + iarr, n);
}
/*----------------------------------------------------------------------------*/

int rd_gse(const RDOM *dom, int iarr, IPNT gs, IPNT ge)
/*<*
 * Get the start and end indices of a specified working (computational virtual) arry in a domain
 *
 * Refer to \ref ra_gse.
 >*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_gse(dom->_s + iarr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int rd_ndim(const RDOM *dom, int iarr, int *ndim)
/*<*
 * Get number of dimensions of a specified array.
 * 
 * @param [in] dom - (const RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [out] ndim - (int *) number of dimensions
 * 
 * @return 0 on successful completion, 
 * if iarr >= narr, return \ref E_BADARRINDEX as in base/include/utils.h.
 *
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_ndim(dom->_s + iarr, ndim);
}
/*----------------------------------------------------------------------------*/

int rd_empty(RDOM *dom, int iarr, int *empty)
/*<*
 * empty query for a specified array in a domain.
 *
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [out] empty - (int *) 0: nonempty, 1: empty
 * @return 0
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_empty(dom->_s + iarr, empty);
}
/*----------------------------------------------------------------------------*/

int rd_setempty(RDOM *dom, int iarr)
/*<*
 * Set a specified working (computational virtual) array in a domain empty.
 *
 * Error arises if iarr >= narr.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setempty(dom->_s + iarr);
}
/*----------------------------------------------------------------------------*/

#ifdef IWAVE_USE_MPI
/*^*/


int rd_setexchangeinfo(RDOM *dom, int iarr, EXCHANGEINFO *einfo)
/*<*
 * Populates exchange info for a specified array in a domain. 
 *
 * Creates MPI_Datatype inside - do not forget to destroy.
 * Refers to \ref ra_setexchangeinfo and \ref IMODEL::ld_r and IMODEL::ld_s.
 * @param[in] dom - (RDOM *) domain pointer
 * @param[in] iarr - (int) array index
 * @param[out] einfo - (EXCHANGEINFO *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setexchangeinfo(dom->_s + iarr, einfo);
}

#endif
/*^*/


/*----------------------------------------------------------------------------*/

int rd_overlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2, int *overlap)
/*<*
 * Checks if two specified working (computational virtual) arrays in two domains overlap.
 * 
 * Refers to \ref ra_overlap.
 * @param [in] dom1, dom2 - (RARR *) domain pointers
 * @param [in] iarr1, iarr2 - (int) array indices 
 * @param [out] overlap - (int *) 0: not overlap, 1: overlap
 * @return 0 on successful completion, else error code as in base/include/utils.h.
>*/
{
    if ( (unsigned int)iarr1 >= (unsigned int)(dom1->narr) ) return E_BADARRINDEX;
    if ( (unsigned int)iarr2 >= (unsigned int)(dom2->narr) ) return E_BADARRINDEX;
	
	return ra_overlap(dom1->_s + iarr1, dom2->_s + iarr2, overlap);
}
/*----------------------------------------------------------------------------*/

int rd_setoverlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2)
/*<*
 * Set the first working (computational virtual) array's dimension info in dom1 to be that of the overlap part 
 * of the two working (computational virtual) arrays in dom1 and dom2.
 *
 * Refer to \ref ra_setoverlap.
 * @param [in,out] dom1 - (RDOM *) domain pointer
 * @param [in] dom2 - (RDOM *) domain pointer
 * @param [in] iarr1, iarr2 - (int) array indices
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 * One of the errors arises when iarr1 >= narr1 or iarr2 >= narr2.
>*/
{
    if ( (unsigned int)iarr1 >= (unsigned int)(dom1->narr) ) return E_BADARRINDEX;
    if ( (unsigned int)iarr2 >= (unsigned int)(dom2->narr) ) return E_BADARRINDEX;
	
	return ra_setoverlap(dom1->_s + iarr1, dom2->_s + iarr2);
}
/*----------------------------------------------------------------------------*/

int rd_a_narr(const RDOM *dom, int *narr)
/*< Get number of arrays in a domain. >*/
{
	*narr = dom->narr;
    return 0;
}
/*----------------------------------------------------------------------------*/
