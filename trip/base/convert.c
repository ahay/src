/* Functions to convert string to a value.

Unlike strtod (and others):
  1) Skips everything (not only white spaces), until finds number start.
  2) Does not change _p, if no conversion can be performed.
  3) In over/underflow cases conversion is performed and corresponding flag.
  
Examples:

strz2double of "abc!2.34Qqqq" will return 2.34 and new string "Qqqq"
--------------------------------------------------------------------------------

strz2int of "abc!2.34Qqqq" will return 2 and new string ".34Qqqq".
Subsequent strz2int will return 34.
*/
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
convert.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
/*^*/

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
/*< strz2char >*/
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
/*< strz2short >*/
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
/*< strz2int >*/
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
/*< strz2long >*/
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
/*< strz2ushort >*/
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
/*< strz2uint >*/
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
/*< strz2ulong >*/
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
/*< strz2double >*/
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
/*< strz2float >*/
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
