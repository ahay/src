/* 
stencil.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "stencil.h"
/*----------------------------------------------------------------------------*/

int mask_setnull(STENCIL_MASK *mask)
{
    memset((void*)mask, 0, sizeof(STENCIL_MASK));

    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_create(STENCIL_MASK *mask, int ip, int ir, int n)
{
    mask_setnull(mask);                 /* empty mask */
    if ( n < 0 ) return E_BADINPUT;

    if ( n > 0 )                        /* allocate memory */
    {
        mask->s = (IPNT*)usermalloc_(n * sizeof(IPNT));
        if ( mask->s == NULL ) return E_ALLOC;
    }
    
    mask->n = n;
    mask->ip = ip;
    mask->ir = ir;

    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_destroy(STENCIL_MASK *mask)
{
    userfree_(mask->s);
    mask_setnull(mask);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_set(STENCIL_MASK *mask, int i, const IPNT ind)
{
    if ( (unsigned int)i >= (unsigned int)(mask->n) ) return E_BADINDEX;
    IASN(mask->s[i], ind);
        
    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_get(STENCIL_MASK *mask, int i, IPNT ind)
{
    if ( (unsigned int)i >= (unsigned int)(mask->n) ) return E_BADINDEX;
    IASN(ind, mask->s[i]);
        
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_setnull(STENCIL *sten)
{
    memset((void*)sten, 0, sizeof(STENCIL));
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_create(STENCIL *sten, int nmask)
{
    int m;
    
    sten_setnull(sten);                 /* empty stencil */
    if ( nmask < 0 ) return E_BADINPUT;

    if ( nmask > 0 )                    /* allocate memory */
    {
        sten->masks = (STENCIL_MASK*)usermalloc_(nmask* sizeof(STENCIL_MASK));
        if ( sten->masks == NULL ) return E_ALLOC;
        for ( m = 0; m < nmask; ++m ) mask_setnull(sten->masks + m);
    }
    
    sten->nmask = nmask;
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_destroy(STENCIL *sten)
{
    int m;
    
    for ( m = 0; m < sten->nmask; ++m ) mask_destroy(sten->masks + m);
    userfree_(sten->masks);
    sten_setnull(sten);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_set(STENCIL *sten, int imask, STENCIL_MASK *mask)
{
    if ( (unsigned int)imask >= (unsigned int)(sten->nmask) ) return E_BADINDEX;
        
    sten->masks[imask] = *mask;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_get(STENCIL *sten, int imask, STENCIL_MASK *mask)
{
    if ( (unsigned int)imask >= (unsigned int)(sten->nmask) ) return E_BADINDEX;
        
    *mask = sten->masks[imask];
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_out(STENCIL *sten, FILE* stream, const char* ind2str_fun(int))
{
    int err, m, i, d;
    STENCIL_MASK mask;
    IPNT ind;
    
    for ( m = 0; m < sten->nmask; ++m )
    {
        err = sten_get(sten, m, &mask);
        if ( err ) return err;
		if ( ind2str_fun )
            printf("mask %d:  ip = %s  ir = %s  n = %d\n", m,
			       ind2str_fun(mask.ip), ind2str_fun(mask.ir), mask.n);
		else
            printf("mask %d:  ip = %d  ir = %d  n = %d\n", m,
				   mask.ip, mask.ir, mask.n);
    
        for ( i = 0; i < mask.n; ++i )
        {
            err = mask_get(&mask, i, ind);
            if ( err ) return err;
            printf(" [");
            for ( d = 0; d < RARR_MAX_NDIM; ++d ) printf(" %d", ind[d]);
            printf(" ]");
        }
		printf("\n");
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
