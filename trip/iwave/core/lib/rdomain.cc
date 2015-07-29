/* 
rdomain.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "rdomain.h"
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
    /*    fprintf(stderr,"rd_a_create_: narr=%d\n",narr); */

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
{
    memset((void*)dom, 0, sizeof(RDOM));    

    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_setnull(RDOM *dom, int iarr)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setnull(dom->_s + iarr);
}
/*----------------------------------------------------------------------------*/

int rd_a_create_s(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])
{
    return rd_a_create_(dom, narr, ndim, dgs, dn, ra_create_s);
}
/*----------------------------------------------------------------------------*/

int rd_a_create_e(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])
{
    return rd_a_create_(dom, narr, ndim, dge, dn, ra_create_e);
}
/*----------------------------------------------------------------------------*/

int rd_a_create(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])
{
    /*  fprintf(stderr,"rd_a_create\n"); */
    return rd_a_create_(dom, narr, ndim, dgs, dge, ra_create);
}
/*----------------------------------------------------------------------------*/

int rd_create_s(RDOM *dom, int ndim, IPNT gs, IPNT n)
{
    return rd_create_(dom, ndim, gs, n, ra_create_s);
}
/*----------------------------------------------------------------------------*/

int rd_create_e(RDOM *dom, int ndim, IPNT ge, IPNT n)
{
    return rd_create_(dom, ndim, ge, n, ra_create_e);
}
/*----------------------------------------------------------------------------*/

int rd_create(RDOM *dom, int ndim, IPNT gs, IPNT ge)
{
    return rd_create_(dom, ndim, gs, ge, ra_create);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare_s(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])
{
    return rd_a_create_(dom, narr, ndim, dgs, dn, ra_declare_s);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare_e(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])
{
    return rd_a_create_(dom, narr, ndim, dge, dn, ra_declare_e);
}
/*----------------------------------------------------------------------------*/

int rd_a_declare(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])
{
    return rd_a_create_(dom, narr, ndim, dgs, dge, ra_declare);
}
/*----------------------------------------------------------------------------*/

int rd_declare_s(RDOM *dom, int ndim, IPNT gs, IPNT n)
{
    return rd_create_(dom, ndim, gs, n, ra_declare_s);
}
/*----------------------------------------------------------------------------*/

int rd_declare_e(RDOM *dom, int ndim, IPNT ge, IPNT n)
{
    return rd_create_(dom, ndim, ge, n, ra_declare_e);
}
/*----------------------------------------------------------------------------*/

int rd_declare(RDOM *dom, int ndim, IPNT gs, IPNT ge)
{
    return rd_create_(dom, ndim, gs, ge, ra_declare);
}
/*----------------------------------------------------------------------------*/
/*
Allocates all arrays.
*/
int rd_a_allocate(RDOM *dom)
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
{
    int err;
    
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    err = ra_allocate(dom->_s + iarr);
    if ( err ) return err;

    return 0;
}
/*----------------------------------------------------------------------------*/

int rd_a_destroy(RDOM *dom)
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
{
    return rd_set_(dom, iarr, gs, n, ra_greset_s);
}
/*----------------------------------------------------------------------------*/

int rd_greset_e(RDOM *dom, int iarr, IPNT ge, IPNT n)
{
    return rd_set_(dom, iarr, ge, n, ra_greset_e);
}
/*----------------------------------------------------------------------------*/

int rd_greset(RDOM *dom, int iarr, IPNT gs, IPNT ge)
{
    return rd_set_(dom, iarr, gs, ge, ra_greset);
}
/*----------------------------------------------------------------------------*/

int rd_offset_s(RDOM *dom, int iarr, IPNT os, IPNT n)
{
    return rd_set_(dom, iarr, os, n, ra_offset_s);
}
/*----------------------------------------------------------------------------*/

int rd_offset_e(RDOM *dom, int iarr, IPNT oe, IPNT n)
{
    return rd_set_(dom, iarr, oe, n, ra_offset_e);
}
/*----------------------------------------------------------------------------*/

int rd_offset(RDOM *dom, int iarr, IPNT os, IPNT oe)
{
    return rd_set_(dom, iarr, os, oe, ra_offset);
}
/*----------------------------------------------------------------------------*/

int rd_a_greset(RDOM *dom, IPNT dgs[], IPNT dge[])
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
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_dump(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_a_print(RDOM *dom, FILE* stream)
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
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_print(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_fprint(RDOM *dom, int iarr, const char *path)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fprint(dom->_s + iarr, path);
}
/*----------------------------------------------------------------------------*/

int rd_write(RDOM *dom, int iarr, FILE* stream)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_write(dom->_s + iarr, stream);
}
/*----------------------------------------------------------------------------*/

int rd_fwrite(RDOM *dom, int iarr, const char *path)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fwrite(dom->_s + iarr, path);
}
/*----------------------------------------------------------------------------*/

int rd_printslice(RDOM *dom, int iarr, FILE* stream, int idim, int li)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_printslice(dom->_s + iarr, stream, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_fprintslice(RDOM *dom, int iarr, const char *path, int idim, int li)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fprintslice(dom->_s + iarr, path, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_writeslice(RDOM *dom, int iarr, FILE* stream, int idim, int li)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_writeslice(dom->_s + iarr, stream, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_fwriteslice(RDOM *dom, int iarr, const char *path, int idim, int li)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_fwriteslice(dom->_s + iarr, path, idim, li);
}
/*----------------------------------------------------------------------------*/

int rd_a_fsprint(RDOM *dom, const char *path)
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
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_get");
    #endif
    
    return ra_get(dom->_s + iarr, li);
}
/*----------------------------------------------------------------------------*/

ireal rd_gget(const RDOM *dom, int iarr, IPNT gi)
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_gget");
    #endif
    
    return ra_gget(dom->_s + iarr, gi);
}
/*----------------------------------------------------------------------------*/

void rd_set(RDOM *dom, int iarr, IPNT li, ireal r)
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_set");
    #endif
    
    ra_set(dom->_s + iarr, li, r);
}
/*----------------------------------------------------------------------------*/

void rd_gset(RDOM *dom, int iarr, IPNT gi, ireal r)
{
    #ifdef CHECK_BOUNDS
    rd_bai(dom, iarr, "rd_gset");
    #endif
    
    ra_gset(dom->_s + iarr, gi, r);
}
/*----------------------------------------------------------------------------*/

int rd_size(RDOM *dom, int iarr, IPNT n)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_size(dom->_s + iarr, n);
}
/*----------------------------------------------------------------------------*/

int rd_a_size(RDOM *dom, int iarr, IPNT n)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_a_size(dom->_s + iarr, n);
}
/*----------------------------------------------------------------------------*/

int rd_gse(const RDOM *dom, int iarr, IPNT gs, IPNT ge)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_gse(dom->_s + iarr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int rd_a_gse(const RDOM *dom, int iarr, IPNT gs, IPNT ge)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_a_gse(dom->_s + iarr, gs, ge);
}
/*----------------------------------------------------------------------------*/

int rd_ndim(const RDOM *dom, int iarr, int *ndim)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_ndim(dom->_s + iarr, ndim);
}
/*----------------------------------------------------------------------------*/

int rd_empty(RDOM *dom, int iarr, int *empty)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_empty(dom->_s + iarr, empty);
}
/*----------------------------------------------------------------------------*/

int rd_setempty(RDOM *dom, int iarr)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setempty(dom->_s + iarr);
}
/*----------------------------------------------------------------------------*/

int rd_setexchangeinfo(RDOM *dom, int iarr, EXCHANGEINFO *einfo)
{
    if ( (unsigned int)iarr >= (unsigned int)(dom->narr) ) return E_BADARRINDEX;
    
    return ra_setexchangeinfo(dom->_s + iarr, einfo);
}
/*----------------------------------------------------------------------------*/

int rd_overlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2, int *overlap)
{
    if ( (unsigned int)iarr1 >= (unsigned int)(dom1->narr) ) return E_BADARRINDEX;
    if ( (unsigned int)iarr2 >= (unsigned int)(dom2->narr) ) return E_BADARRINDEX;
	
	return ra_overlap(dom1->_s + iarr1, dom2->_s + iarr2, overlap);
}
/*----------------------------------------------------------------------------*/

int rd_setoverlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2)
{
    if ( (unsigned int)iarr1 >= (unsigned int)(dom1->narr) ) return E_BADARRINDEX;
    if ( (unsigned int)iarr2 >= (unsigned int)(dom2->narr) ) return E_BADARRINDEX;
	
	return ra_setoverlap(dom1->_s + iarr1, dom2->_s + iarr2);
}
/*----------------------------------------------------------------------------*/

int rd_a_narr(const RDOM *dom, int *narr)
{
	*narr = dom->narr;
    return 0;
}

/*----------------------------------------------------------------------------*/

int rd_a_inner(RDOM const * dom1, RDOM const * dom2, ireal * ip) {
  int a;
  ireal ipp = REAL_ZERO;
  *ip = REAL_ZERO;
  int err = 0;
  
  if (dom1->narr != dom2->narr) {
    err=E_BADINPUT;
    return err;
  }
  for ( a = 0; a < dom1->narr; ++a ) {
    if (!(err = ra_a_inner(dom1->_s + a, dom2->_s + a, &ipp))) 
      *ip += ipp;
    else return err;
  }
  
  return err;
}

/*----------------------------------------------------------------------------*/
int rd_a_zero(RDOM * dom) {

    int a;
    int err = 0;
    
    for ( a = 0; a < dom->narr; ++a )
    {
      err = ra_a_zero(dom->_s + a);
      if (err) return err;
    }
    
    return err;
}

/*----------------------------------------------------------------------------*/

int rd_a_scale(RDOM * dom, int iarr, ireal fac) {
  int err = 0;  
  err = ra_a_scale(dom->_s + iarr, fac); 
  return err;
}

