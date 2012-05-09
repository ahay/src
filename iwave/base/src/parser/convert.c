/* 
convert.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "convert.h"
/*----------------------------------------------------------------------------*/
/*
Checks for an floating number start.
Returns true in the following cases: 
    [d], [+d], [-d], [.d], [+.d], [-.d].
Where d is a digit.
Returns false otherwise.

char *str :  null-terminated string.
*/
static inline int flt_start(char *str)
{
    /* TODO: check nan here */
    if ( (*str == '+') || (*str == '-') ) ++str;
    /* TODO: check inf here */
    if ( *str == '.' ) ++str;
    return isdigit(*str);
}
/*----------------------------------------------------------------------------*/
/*
Checks for a float number start.
Returns true in the following cases: [d], [+d], [-d]. Where d is a digit.
Returns false otherwise.

char *str :  null-terminated string.
*/
static inline int int_start(char *str)
{
    if ( (*str == '+') || (*str == '-') ) ++str;
    return isdigit(*str);
}
/*----------------------------------------------------------------------------*/

int strz2char(SIZEDSTRING str, char *p, SIZEDSTRING *strend)
{
    if ( str.n == 0L ) 
    {
        if ( strend != NULL ) *strend = str;
        return E_PARSECONVERT;
    }
    if ( p != NULL ) *p = *(str.s);
    if ( strend != NULL )
    {
        strend->n = str.n - 1L;
        strend->s = str.s + 1;
    }
    return 0;
}
/*----------------------------------------------------------------------------*/

int strz2short(SIZEDSTRING str, short *p, SIZEDSTRING *strend)
{
    long l;
    int err;

    err = strz2long(str, &l, strend);          /* convert to long */
    if ( err && (err != E_RANGE) ) return err; /* conversion error */

    if ( l > (long)SHRT_MAX )       /* overflow */
    {
        if ( p != NULL  ) *p = SHRT_MAX;
        err = E_RANGE;
    }
    else if ( l < (long)SHRT_MIN )  /* underflow */
    {
        if ( p != NULL  ) *p = SHRT_MIN;
        err = E_RANGE;
    }
    else if ( p != NULL ) *p = (short)l;  /* no error */
    return err;
}
/*----------------------------------------------------------------------------*/

int strz2int(SIZEDSTRING str, int *p, SIZEDSTRING *strend)
{
    long l;
    int err;

    err = strz2long(str, &l, strend);          /* convert to long */
    if ( err && (err != E_RANGE) ) return err; /* conversion error */

    if ( l > (long)INT_MAX )       /* overflow */
    {
        if ( p != NULL  ) *p = INT_MAX;
        err = E_RANGE;
    }
    else if ( l < (long)INT_MIN )  /* underflow */
    {
        if ( p != NULL  ) *p = INT_MIN;
        err = E_RANGE;
    }
    else if ( p != NULL ) *p = (int)l;  /* no error */
    return err;
}
/*----------------------------------------------------------------------------*/

int strz2long(SIZEDSTRING str, long *p, SIZEDSTRING *strend)
{
    long l;
    SIZEDSTRING end;
    if ( strend == NULL ) strend = &end;
    
    /* skip symbols */
    while ( str.n > 0L )
    {
        if ( int_start(str.s) ) break;
        ++(str.s);
        --(str.n);
    }
    
    /* convert */
    errno = 0;
    l = strtol(str.s, &(strend->s), 10);
    strend->n = str.n - (strend->s - str.s);
    
    /* check error  */
    if ( str.s == strend->s ) return E_PARSECONVERT; /* cannot convert */
    if ( p != NULL ) *p = l;
    if ( errno == ERANGE ) return E_RANGE;
    return 0;
}
/*----------------------------------------------------------------------------*/

int strz2ushort(SIZEDSTRING str, unsigned short *p, SIZEDSTRING *strend)
{
    unsigned long l;
    int err;

    err = strz2ulong(str, &l, strend);           /* convert to long */
    if ( err && (err != E_RANGE) ) return err;   /* conversion error */

    if ( l > (unsigned long)USHRT_MAX )          /* overflow */
    {
        if ( p != NULL  ) *p = USHRT_MAX;
        err = E_RANGE;
    }
    else if ( p != NULL ) *p = (unsigned short)l;/* no error */
    return err;
}
/*----------------------------------------------------------------------------*/

int strz2uint(SIZEDSTRING str, unsigned int *p, SIZEDSTRING *strend)
{
    unsigned long l;
    int err;

    err = strz2ulong(str, &l, strend);          /* convert to long */
    if ( err && (err != E_RANGE) ) return err;  /* conversion error */

    if ( l > (unsigned long)UINT_MAX )          /* overflow */
    {
        if ( p != NULL  ) *p = UINT_MAX;
        err = E_RANGE;
    }
    else if ( p != NULL ) *p = (unsigned int)l; /* no error */
    return err;
}
/*----------------------------------------------------------------------------*/

int strz2ulong(SIZEDSTRING str, unsigned long *p, SIZEDSTRING *strend)
{
    unsigned long l;
    SIZEDSTRING end;
    if ( strend == NULL ) strend = &end;
    
    /* skip symbols */
    while ( str.n > 0L )
    {
        if ( int_start(str.s) ) break;
        ++(str.s);
        --(str.n);
    }
    
    /* convert */
    errno = 0;
    l = strtoul(str.s, &(strend->s), 10);
    strend->n = str.n - (strend->s - str.s);
    
    /* check error  */
    if ( str.s == strend->s ) return E_PARSECONVERT; /* cannot convert */
    if ( p != NULL ) *p = l;
    if ( errno == ERANGE ) return E_RANGE;
    return 0;
}
/*----------------------------------------------------------------------------*/
int strz2double(SIZEDSTRING str, double *p, SIZEDSTRING *strend)
{
    double d;
    SIZEDSTRING end;
    if ( strend == NULL ) strend = &end;
    
    /* skip symbols */
    while ( str.n > 0L )
    {
        if ( flt_start(str.s) ) break;
        ++(str.s);
        --(str.n);
    }
    
    /* convert */
    errno = 0;
    d = strtod(str.s, &(strend->s));
    strend->n = str.n - (strend->s - str.s);
    
    /* check error  */
    if ( str.s == strend->s ) return E_PARSECONVERT; /* cannot convert */
    if ( p != NULL ) *p = d;
    if ( errno == ERANGE ) return (d == 0.0) ? E_UNDERFLOW : E_OVERFLOW;
    return 0;
}
/*----------------------------------------------------------------------------*/

int strz2float(SIZEDSTRING str, float *p, SIZEDSTRING *strend)
{
#if __STDC_VERSION__ >= 199901L
    float f;
    SIZEDSTRING end;
    if ( strend == NULL ) strend = &end;
    
    /* skip symbols */
    while ( str.n > 0L )
    {
        if ( flt_start(str.s) ) break;
        ++(str.s);
        --(str.n);
    }
    
    /* convert */
    errno = 0;
    f = strtof(str.s, &(strend->s));
    strend->n = str.n - (strend->s - str.s);
    
    /* check error  */
    if ( str.s == strend->s ) return E_PARSECONVERT; /* cannot convert */
    if ( p != NULL ) *p = f;
    if ( errno == ERANGE ) return (f == 0.0f) ? E_UNDERFLOW : E_OVERFLOW;
    return 0;
#else
    double d;
    int err;
    
    err = strz2double(str, &d, strend);

    if ( err && (err != E_UNDERFLOW) && (err != E_OVERFLOW) ) return err;
        
    if ( d > (double)FLT_MAX )          /* overflow */
    {
        if ( p != NULL  ) *p = FLT_MAX;
        err = E_OVERFLOW;
    }
    else if ( d < (double)(-FLT_MAX) )  /* overflow */
    {
        if ( p != NULL  ) *p = -FLT_MAX;
        err = E_OVERFLOW;
    }
    else
    {
        if ( p != NULL  ) *p = (float)d;
        if ( (d != 0.0) && (((float)d) == 0.0f) ) err = E_UNDERFLOW;
    }
    return err;
#endif
}
/*----------------------------------------------------------------------------*/
