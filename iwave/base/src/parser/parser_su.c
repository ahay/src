/* 
parser_su.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "parser.h"
#include "parser_su.h"
/*----------------------------------------------------------------------------*/

static PARARRAY m_parr = {0, NULL, NULL, NULL};
/*----------------------------------------------------------------------------*/
/*
Converts char type to int type.
*/
static int ps_type2int(char *type, int *itype, size_t *size)
{
    switch ( *type )
    {
        case 'i':
            *itype = DT_INT;
            *size = sizeof(int);
            break;
        case 'p':
            *itype = DT_UINT;
            *size = sizeof(unsigned int);
            break;
        case 'h':
            *itype = DT_SHORT;
            *size = sizeof(short);
            break;
        case 'u':
            *itype = DT_USHORT;
            *size = sizeof(unsigned short);
            break;
        case 'l':
            *itype = DT_LONG;
            *size = sizeof(long);
            break;
        case 'v':
            *itype = DT_ULONG;
            *size = sizeof(unsigned long);
            break;
        case 'f':
            *itype = DT_FLOAT;
            *size = sizeof(float);
            break;
        case 'd':
            *itype = DT_DOUBLE;
            *size = sizeof(double);
            break;
        case 'r':
            *itype = DT_REAL;
            *size = sizeof(ireal);
            break;
        case 's':
            *itype = DT_CSTRING;
            *size = 0;
            break;
        default:
            return E_BADINPUT;
    }
    return 0;
}
/*----------------------------------------------------------------------------*/

void ps_initargs(int argc, char **argv)
{
    PARARRAY parr_arg;                  /* array of arguments */
    SIZEDSTRING parfile;                /* name of parameter file */
    
    ps_destroy(&m_parr);
    
    /* create paramater list from arguments (skip program name) */
    if ( ps_createargs(&parr_arg, argc - 1, argv + 1) ) return;
    
    /* if parameter file missing, use arguments */    
    if ( ps_getval(parr_arg, "par", 0, &parfile) )
    {
        m_parr = parr_arg;
        return;
    }
    
    /* load from given parameter file */
    ps_createfile(&m_parr, parfile.s);
    ps_destroy(&parr_arg);
}
/*----------------------------------------------------------------------------*/

void ps_destroyargs()
{
    ps_destroy(&m_parr);
}
/*----------------------------------------------------------------------------*/

int ps_countparname(char *name)
{
    return ps_countname(m_parr, name);
}
/*----------------------------------------------------------------------------*/

int ps_countparval(char *name)
{
    return ps_countnparval(0, name);
}
/*----------------------------------------------------------------------------*/

int ps_countnparval(int n, char *name)
{
    SIZEDSTRING val;                    /* parameter value */
    int count;
    long i;
    
    /* n starts from, 0 means last element */
    if ( n == 0 ) n = ps_countname(m_parr, name);
    
    /* get parameter value */
    if ( ps_getval(m_parr, name, n - 1, &val) ) return 0;
        
    /* count ',' symbols */
    count = 1;
    for ( i = 0L; i < val.n; ++i ) if ( val.s[i] == ',' ) ++count;
        
    return count;
}
/*----------------------------------------------------------------------------*/

void* ps_getpararray()
{
    return (void*)(&m_parr);
}
/*----------------------------------------------------------------------------*/

int ps_getparint(char *name, int *p)
{
    return ps_getnparint(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparuint(char *name, unsigned int *p)
{
    return ps_getnparuint(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparshort(char *name, short *p)
{
    return ps_getnparshort(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparushort(char *name, unsigned short *p)
{
    return ps_getnparushort(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparlong(char *name, long *p)
{
    return ps_getnparlong(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparulong(char *name, unsigned long *p)
{
    return ps_getnparulong(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparfloat(char *name, float *p)
{
    return ps_getnparfloat(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getpardouble(char *name, double *p)
{
    return ps_getnpardouble(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getparreal(char *name, ireal *p)
{
    return ps_getnparreal(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparint(int n, char *name, int *p)
{
    return ps_getnpar(n, name, "i", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparuint(int n, char *name, unsigned int *p)
{
    return ps_getnpar(n, name, "p", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparshort(int n, char *name, short *p)
{
    return ps_getnpar(n, name, "h", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparushort(int n, char *name, unsigned short *p)
{
    return ps_getnpar(n, name, "u", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparlong(int n, char *name, long *p)
{
    return ps_getnpar(n, name, "l", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparulong(int n, char *name, unsigned long *p)
{
    return ps_getnpar(n, name, "v", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparfloat(int n, char *name, float *p)
{
    return ps_getnpar(n, name, "f", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnpardouble(int n, char *name, double *p)
{
    return ps_getnpar(n, name, "d", p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparreal(int n, char *name, ireal *p)
{
    return ps_getnpar(n, name, "r", p);
}
/*----------------------------------------------------------------------------*/

void ps_getPar(char *name, char *type, void *p)
{
    ps_getnpar(0, name, type, p);
}
/*----------------------------------------------------------------------------*/

int ps_getnpar(int n, char *name, char *type, void *p)
{
    int maxcount, i, itype;
    size_t shift;
    SIZEDSTRING val;                    /* parameter value */
    
    /* number of values */
    maxcount = ps_countnparval(n, name);
    if ( maxcount <= 0 ) return 0;
    
    /* n starts from, 0 means last element */
    if ( n == 0 ) n = ps_countname(m_parr, name);
    
    /* get parameter value */
    if ( ps_getval(m_parr, name, n - 1, &val) ) return 0;

    if ( *type == 's' ) maxcount = 1;   /* string array processed separately */
        
    /* get DT_TYPE and type size */
    if ( ps_type2int(type, &itype, &shift) ) return 0;
        
    for ( i = 0; i < maxcount; ++i )
    {
        if ( ps_val2type(val, itype, p, &val) ) break;
        p = (void*)((char*)p + shift);
    }
        
    return i;
}
/*----------------------------------------------------------------------------*/

int ps_getparstring(char *name, char **p)
{
    return ps_getnparstring(0, name, p);
}
/*----------------------------------------------------------------------------*/

int ps_getnparstring(int n, char *name, char **p)
{
    return ps_getnpar(n, name, "s", p);
}
/*----------------------------------------------------------------------------*/
