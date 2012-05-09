/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/* HDRPKGE: $Revision: 1.11 $ ; $Date: 1996/09/09 17:13:34 $	*/

/*********************** self documentation **********************/
/************************************************************************** 
HDRPKGE - routines to access the SEGY header via the hdr structure.
 
gethval		get a trace header word by index
puthval		put a trace header word by index
getbhval	get a binary header word by index
putbhval	put a binary header word by index
gethdval	get a trace header word by name
puthdval	put a trace header word by name


hdtype		get the data type of a trace header word by name
getkey		get the name of a trace header word from its index

getindex	get the index of a trace header word from the name

swaphval	swap the trace header words by index
swapbhval	swap the binary header words by index
gettapehval	get a tape trace header word by index
puttapehval	put a tape trace header word by index
gettapebhval	get a tape binary header word by index
puttapebhval	put a tape binary header word by index
printheader	display non-null header field values

*************************************************************************** 
Function Prototypes:
void gethval(const segy *tr, int index, Value *valp);
void puthval(segy *tr, int index, Value *valp);
void putbhval(bhed *bh, int index, Value *valp);
void getbhval(const bhed *bh, int index, Value *valp);
void gethdval(const segy *tr, char *key, Value *valp);
void puthdval(segy *tr, char *key, Value *valp);
char *hdtype(const char *key);
char *getkey(const int index);
int getindex(const char *key);
void swaphval(segy *tr, int index);
void swapbhval(bhed *bh, int index);
void gettapehval(tapesegy *tapetr, int index, Value *valp);
void puttapehval(tapesegy *tapetr, int index, Value *valp);
void gettapebhval(tapebhed *tapetr, int index, Value *valp);
void puttapebhval(tapebhed *tapetr, int index, Value *valp);
void printheader(const segy *tp);

*************************************************************************** 
Notes:
This package includes only those routines that directly access
the "hdr" or "bhdr" structures.  It does not include routines
such as printfval, printftype, printfhead that use the routines
in this package to indirectly access these structures.

Note that while gethdval and puthdval are more convenient to use
than gethval and puthval, they incur an inefficiency in the
common case of iterating code over a set of traces with a fixed
key or keys.  In such cases, it is advisable to set the index
or indices outside the loop using getindex.

swaphval:
Byte-swapping is needed for converting SU data from big-endian to little-
endian formats, and vice versa. The swap_.... subroutines are based
on subroutines provided by Jens Hartmann of the Institut fur Geophysik
in Hamburg. These are found in .../cwp/lib/swapbyte.c.
**************************************************************************
Authors: SEP: Einar Kjartansson	CWP: Jack Cohen, Shuki Ronen
swaphval: CWP: John Stockwell
*************************************************************************/
/**************** end self doc ********************************/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "hdr.h"
#include "bhdr.h"
#include "tapesegy.h"
#include "tapehdr.h"
#include "tapebhdr.h"

void gethval(const segy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(valp->s, tp + hdr[index].offs);  break;
	case 'h': valp->h = *((short*)  (tp + hdr[index].offs)); break;
	case 'u': valp->u = *((unsigned short*) (tp + hdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (tp + hdr[index].offs)); break;
	case 'p': valp->p = *((unsigned int*)   (tp + hdr[index].offs)); break;
	case 'l': valp->l = *((long*)   (tp + hdr[index].offs)); break;
	case 'v': valp->v = *((unsigned long*)  (tp + hdr[index].offs)); break;
	case 'f': valp->f = *((float*)  (tp + hdr[index].offs)); break;
	case 'd': valp->d = *((double*) (tp + hdr[index].offs)); break;
	default: suerr("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void puthval(segy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(tp + hdr[index].offs, valp->s);  break;
	case 'h': *((short*)  (tp + hdr[index].offs)) = valp->h; break;
	case 'u': *((unsigned short*) (tp + hdr[index].offs)) = valp->u; break;
	case 'i': *((int*)   (tp + hdr[index].offs)) = valp->i; break;
	case 'p': *((unsigned int*)   (tp + hdr[index].offs)) = valp->p; break;
	case 'l': *((long*)   (tp + hdr[index].offs)) = valp->l; break;
	case 'v': *((unsigned long*)  (tp + hdr[index].offs)) = valp->v; break;
	case 'f': *((float*)  (tp + hdr[index].offs)) = valp->f; break;
	case 'd': *((double*) (tp + hdr[index].offs)) = valp->d; break;
	default: suerr("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void getbhval(const bhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(bhdr[index].type)) {
	case 'h': valp->h = *((short*) (bhp + bhdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (bhp + bhdr[index].offs)); break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void putbhval(bhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(bhdr[index].type)) {
	case 'h': *((short*) (bhp + bhdr[index].offs)) = valp->h; break;
	case 'i': *((int*)   (bhp + bhdr[index].offs)) = valp->i; break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void gethdval(const segy *tr, char *key, Value *valp)
{
	int index = getindex(key);
	char *tp = (char*) tr;

	if ( -1 == (index))
		suerr("%s: key word not in segy.h: '%s'", __FILE__, key);

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(valp->s, tp + hdr[index].offs);  break;
	case 'h': valp->h = *((short*)  (tp + hdr[index].offs)); break;
	case 'u': valp->u = *((unsigned short*) (tp + hdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (tp + hdr[index].offs)); break;
	case 'p': valp->p = *((unsigned int*)   (tp + hdr[index].offs)); break;
	case 'l': valp->l = *((long*)   (tp + hdr[index].offs)); break;
	case 'v': valp->v = *((unsigned long*)  (tp + hdr[index].offs)); break;
	case 'f': valp->f = *((float*)  (tp + hdr[index].offs)); break;
	case 'd': valp->d = *((double*) (tp + hdr[index].offs)); break;
	default: suerr("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void puthdval(segy *tr, char *key, Value *valp)
{
	int index = getindex(key);
	char *tp = (char*) tr;

	if ( -1 == (index))
		suerr("%s: key word not in segy.h: '%s'", __FILE__, key);

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(tp + hdr[index].offs, valp->s);  break;
	case 'h': *((short*)  (tp + hdr[index].offs)) = valp->h; break;
	case 'u': *((unsigned short*) (tp + hdr[index].offs)) = valp->u; break;
	case 'i': *((int*)   (tp + hdr[index].offs)) = valp->i; break;
	case 'p': *((unsigned int*)   (tp + hdr[index].offs)) = valp->p; break;
	case 'l': *((long*)   (tp + hdr[index].offs)) = valp->l; break;
	case 'v': *((unsigned long*)  (tp + hdr[index].offs)) = valp->v; break;
	case 'f': *((float*)  (tp + hdr[index].offs)) = valp->f; break;
	case 'd': *((double*) (tp + hdr[index].offs)) = valp->d; break;
	default: suerr("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


char *hdtype(const char *key)
{
	int index = getindex(key);

	if (-1 == (index))
		suerr("%s: key word not in segy.h: '%s'", __FILE__, key);

	return hdr[index].type;
}


char *getkey(const int index)
{
	return (index < SU_NKEYS && index >= 0) ? hdr[index].key : NULL;
}


int getindex(const char *key)	/* get index for this key */
{
	register int i;

	for (i = 0; i < SU_NKEYS; i++)
		if (STREQ(hdr[i].key, key))
			return i;	/* key found */

	/* not found */
	return -1;
}


void swaphval(segy *tr, int index)
{
	register char *tp= (char *) tr;

        switch(*(hdr[index].type)) {
        case 'h': swap_short_2((short*)(tp + hdr[index].offs));
	break;
        case 'u': swap_u_short_2((unsigned short*)(tp + hdr[index].offs));
	break;
        case 'i': swap_int_4((int*)(tp + hdr[index].offs));
	break;
        case 'p': swap_u_int_4((unsigned int*)(tp + hdr[index].offs));
	break;
        case 'l': swap_long_4((long*)(tp + hdr[index].offs));
	break;
        case 'v': swap_u_long_4((unsigned long*)(tp + hdr[index].offs));
	break;
        case 'f': swap_float_4((float*)(tp + hdr[index].offs));
	break;
        case 'd': swap_double_8((double*)(tp + hdr[index].offs));
	break;
        default: suerr("%s: %s: unsupported data type", __FILE__, __LINE__);
	break;
        }

        return;
}


void swapbhval(bhed *bh, int index)
{
	register char *bhp= (char *) bh;

        switch(*(bhdr[index].type)) {
        case 'h': swap_short_2((short*)(bhp + bhdr[index].offs)); break;
        case 'i': swap_int_4((int*)(bhp + bhdr[index].offs)); break;
        default: suerr("%s: %s: unsupported data type", __FILE__, __LINE__);
	break;
        }

        return;
}


void gettapehval(const tapesegy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapehdr[index].type)) {
	case 'U': valp->h = (short) *((short*) (tp + tapehdr[index].offs));
	break;
	case 'P': valp->i = (int) *((int*) (tp + tapehdr[index].offs));
	break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void puttapehval(tapesegy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapehdr[index].type)) {
	case 'U': *((short*) (tp + tapehdr[index].offs)) = valp->h; break;
	case 'P': *((int*)   (tp + tapehdr[index].offs)) = valp->i; break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void gettapebhval(const tapebhed *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapebhdr[index].type)) {
	case 'U': valp->h = (short) *((short*) (tp + tapebhdr[index].offs));
	break;
	case 'P': valp->i = (int) *((int*) (tp + tapebhdr[index].offs));
	break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void puttapebhval(tapebhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(tapebhdr[index].type)) {
	case 'U': *((short*) (bhp + tapebhdr[index].offs)) = valp->h; break;
	case 'P': *((int*)   (bhp + tapebhdr[index].offs)) = valp->i; break;
	default: suerr("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}



/* Display non-null header field values */
void printheader(const segy *tp)
{
	int i;			/* index over header fields		*/
	int j;			/* index over non-null header fields	*/
	Value val;		/* value in header field		*/
	cwp_String type;	/* ... its data type			*/
	cwp_String key;		/* ... the name of the header field	*/
	Value zeroval;		 /* zero value to compare with		*/

	zeroval.l = 0;
	j = 0;
	for (i = 0; i < SU_NKEYS; i++) {
		gethval(tp, i, &val);
		key = getkey(i);
		type = hdtype(key);
		if (valcmp(type, val, zeroval)) { /* not equal to zero */
			(void) printf(" %s=", key);
			printfval(type, val);
			if ((++j % 6) == 0) putchar('\n');
		}
	}
	putchar('\n');

	return;
}
