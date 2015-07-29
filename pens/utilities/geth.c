/*
 *
 *  source file:   ./filters/loclib/geth.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Martin Karrenbach ifdef cray
 */

#include <stdio.h>

short geth(FILE *iop)
/*< get short >*/
{
    union {
	unsigned short w;
	short s;
    } r;

    r.w = getc(iop);
    r.w += (getc(iop) << 8);
#ifdef CRAY
    if( r.s > 0x8000 ) r.s = r.s - 0x10000;
#endif
    if (feof(iop))
	return(EOF);
    return((int)r.s);
}
