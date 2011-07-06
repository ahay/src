/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
ATOPKGE - convert ascii to arithmetic and with error checking

 
eatoh		ascii to short
eatou		ascii to unsigned short
eatoi		ascii to int
eatop		ascii to unsigned
eatol		ascii to long
eatov		ascii to unsigned long
eatof		ascii to float
eatod		ascii to double

****************************************************************************
Function Prototypes:
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

****************************************************************************
Input:
s		string 

Returned:	type indicated
 
****************************************************************************
Notes:
Each of these routines acts like atoi, but has error checking:

This is a major revision of the tedious code required before
vendors implemented the ANSI C strtol, strtoul and strtod.

In addition to the size checks for each integer type, a
specific test on errno is required.  For example, INT_MAX
may (and probably does) equal LONG_MAX.  In this case,
if fed a number exceeding INT_MAX (and LONG_MAX), strtol
will make a quiet return with the wrong answer and it is up
to the user to check if errno == ERANGE.

Size limits are machine dependent and are read from the
ANSI C include files limits.h and float.h.

Bug Report: With NeXT c and Gnucc, when x > DBL_MAX (x <-DBL_MAX),
the return value from strtod was +Infinity (-Infinity), not HUGE_VAL
and more important, errno was not set to ERANGE.  To cope with this,
I put explicit size checks in eatod (which would not be needed if
errno were set as it should be in ANSI C.    jkc 01/29/94

On IBM RS6000, the return value from strtod was +-Inf on
overflow, but errno was set correctly.

****************************************************************************
References:
For old code:
Plum: Reliable Data Structures in C, p. 2-17.
Kernighan and Ritchie: The C Programming Language, p. 58.

CWP: Jack K. Cohen, Brian Sumner
 
For new code:
ANSI C routines with a little help from Jack

****************************************************************************
Author: Jack Cohen, Center for Wave Phenomena, 1994.
***************************************************************************/
/**************** end self doc ********************************/

#include <stdlib.h>
#include <limits.h>

#include <rsf.h>

#include "atopkge.h"

#ifndef SUN_A
#include <errno.h>

short eatoh(char *s)
/*< eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} >*/
{
    long n = strtol(s, NULL, 10);
	
    if ( (n > SHRT_MAX) || (n < SHRT_MIN) || (errno == ERANGE) )
	sf_error("%s: eatoh: overflow", __FILE__);

    return (short) n;
}

unsigned short eatou(char *s)
/*< eatou - convert string s to unsigned short integer {0:USHRT_MAX} >*/
{
    unsigned long n = strtoul(s, NULL, 10);

    if ( (n > USHRT_MAX) || (errno == ERANGE) )
	sf_error("%s: eatou: overflow", __FILE__);

    return (unsigned short) n;
}

int eatoi(char *s)
/*< eatoi - convert string s to integer {INT_MIN:INT_MAX} >*/
{
    long n = strtol(s, NULL, 10);

    if ( (n > INT_MAX) || (n < INT_MIN) || (errno == ERANGE) )
	sf_error("%s: eatoi: overflow", __FILE__);

    return (int) n;
}

unsigned int eatop(char *s)
/*< eatop - convert string s to unsigned integer {0:UINT_MAX} >*/
{
    unsigned long n = strtoul(s, NULL, 10);

    if ( (n > UINT_MAX) || (errno == ERANGE) )
	sf_error("%s: eatop: overflow", __FILE__);

    return (unsigned int) n;
}

long eatol(char *s)
/*< eatol - convert string s to long integer {LONG_MIN:LONG_MAX} >*/
{
    long n = strtol(s, NULL, 10);

    if (errno == ERANGE)
	sf_error("%s: eatol: overflow", __FILE__);

    return n;
}

unsigned long eatov(char *s)
/*< eatov - convert string s to unsigned long {0:ULONG_MAX} >*/
{
    unsigned long n = strtoul(s, NULL, 10);

    if (errno == ERANGE)
	sf_error("%s: eatov: overflow", __FILE__);

    return n;
}

float eatof(char *s)
/*< eatof - convert string s to float {-FLT_MAX:FLT_MAX} >*/
{
    float x = strtod(s, NULL);

    if ( (x > FLT_MAX) || (x < -FLT_MAX) || (errno == ERANGE) )
	sf_error("%s: eatof: overflow", __FILE__);

    return (float) x;
}

double eatod(char *s)
/*< eatod - convert string s to double {-DBL_MAX:DBL_MAX} >*/
{
    double x = strtod(s, NULL);

    /* errno == ERANGE suffices if compiler sets errno on overflow */
    if ( (errno == ERANGE) || (x > DBL_MAX) || (x < -DBL_MAX) )
	sf_error("%s: eatod: overflow", __FILE__);

    return x;
}

#else   /* if SUN_A is defined */
/* old one resurrected because SUN 2 didn't have strtoul */

/* atopkge - convert ascii to arithmetic and with error checking
 *
 * eatoh	- ascii to short
 * eatou	- ascii to unsigned short
 * eatoi	- ascii to int
 * eatop	- ascii to unsigned
 * eatol	- ascii to long
 * eatov	- ascii to unsigned long
 * eatof	- ascii to float (dummy sub)
 * eatod	- ascii to double (dummy sub)
 *
 * Returns:
 *	eatoh: short int
 *	eatou: unsigned short int
 *	eatoi: int
 *	eatop: unsigned int
 *	eatol: long int
 *	eatov: unsigned long int
 *	eatof: float
 *	eatod: double
 *
 * Synopsis:
 *	short eatoh(s)
 *	char s[];
 *
 *	unsigned short eatou(s)
 *	char s[];
 *
 *	int eatoi(s)
 *	char s[];
 *
 *	unsigned eatop(s)
 *	char s[];
 *
 *	long eatol(s)
 *	char s[];
 *
 *	unsigned long eatov(s)
 *	char s[];
 *
 *	float eatof(s)
 *	char s[];
 *
 *	double eatod(s)
 *	char s[];
 *
 * Notes:
 *	I haven't a clue as to how to write eatof and eatod, but when
 *	vendors come up to snuff on the ANSI C prescribed error returns
 *	for strtod, it'll be a piece of cake.  And when strtol, strtoul
 *	are likewise properly implemented, the remaining routines in this
 *	package will simplify materially.  For now, eatof and eatod are
 *	just place holders that don't really check for errors.
 *
 *	The overflow error check on a type that fills an unsigned long
 *	is different and a bit slower than the others.  Still, it might
 *      be better to use it in eatou and eatop as well and avoid the
 *	(possible) additional function call.
 *
 *	The code relies on the fact that converting unsigned to signed
 *	has no surprises for numbers in the lower half of the range.
 *
 *	Size limits on the integer data types are machine dependent and
 *      are read from the file limits.h.
 *
 * Credits:
 *	Plum: Reliable Data Structures in C, p. 2-17.
 *	Kernighan and Ritchie: The C Programming Language, p. 58.
 *	CWP: Jack, Brian
 *
 *
 */

/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX}    *
 * We store the absolute value of the converted string in an        *
 * unsigned long so we can test it for overflow.                    */
short eatoh(char *s)
{
    unsigned long n;
    int i;
    short sign = 1;
    long eatol();

    for (i = 0; isspace(s[i]); ++i) ;	/* skip white space */

    if (s[i] == '+' || s[i] == '-') {
	sign = (s[i++] == '+') ? 1 : -1;
    }

    for (n = 0; isdigit(s[i]) && n <= SHRT_MAX/10; ++i) {
	n = 10 * n + (s[i] - '0');
    }

    if ((sign ==  1) && (n > SHRT_MAX) ||
	(sign == -1) && (n > SHRT_MIN) || isdigit(s[i]))
	sf_error("%s: eatoh: overflow", __FILE__);

    return  sign * (short) n;
}


/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} *
 * If USHRT_MAX < ULONG_MAX, we can temporarily fit the converted   *
 * number in an unsigned long with room to check for overflow       *
 * condition.  If not, we forward the string to the unsigned long   *
 * routine.                                                         */
unsigned short eatou(char *s)
{
    unsigned long n;
    int i;
    unsigned long eatov();

    if (USHRT_MAX == ULONG_MAX)  return (unsigned short) eatov(s);

    for (i = 0; isspace(s[i]); ++i) ;  /* skip white space */

    if (s[i] == '-')
	sf_error("%s: eatou: saw negative number", __FILE__);

    for (n = 0; isdigit(s[i]) && n <= USHRT_MAX/10; ++i) {
	n = 10 * n + (s[i] - '0');
    }
    if (n > USHRT_MAX || isdigit(s[i]))
	sf_error("%s: eatou: overflow", __FILE__);

    return (unsigned short) n;
}


/* eatoi - convert string s to short integer {INT_MIN:INT_MAX}    *
 * The logic is the same as for eatou with INT_MAX replacing      *
 * SHRT_MAX and INT_MIN replacing SHRT_MIN.                       */
int eatoi(char *s)
{
    unsigned long n;
    int i;
    int sign = 1;
    long eatol();

    if (INT_MAX == LONG_MAX) return (int) eatol(s);

    for (i = 0; isspace(s[i]); ++i) ;	/* skip white space */

    if (s[i] == '+' || s[i] == '-') {
	sign = (s[i++] == '+') ? 1 : -1;
    }

    for (n = 0; isdigit(s[i]) && n <= INT_MAX/10; ++i) {
	n = 10 * n + (s[i] - '0');
    }

    if ((sign ==  1) && (n > INT_MAX) ||
	(sign == -1) && (n > INT_MIN) || isdigit(s[i]))
	sf_error("%s: eatoi: overflow", __FILE__);

    return  sign * (int) n;
}


/* eatop - convert string s to unsigned integer {0:UINT_MAX}        *
 * The logic is the same as for eatou with UINT_MAX replacing       *
 * USHRT_MAX.                                                       */
unsigned int eatop(char *s)
{
    unsigned long n;
    int i;
    unsigned long eatov();

    if (UINT_MAX == ULONG_MAX) return((unsigned int) eatov(s));

    for (i = 0; isspace(s[i]); ++i) ;  /* skip white space */

    if (s[i] == '-')
	sf_error("%s: eatop: saw negative number", __FILE__);

    for (n = 0; isdigit(s[i]) && n <= UINT_MAX/10; ++i) {
	n = 10 * n + (s[i] - '0');
    }
    if (n > UINT_MAX || isdigit(s[i]))
	sf_error("%s: eatop: overflow", __FILE__);

    return (unsigned int) n;
}


/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX}     *
 * We store the absolute value of the converted string in an        *
 * unsigned long so we can test it for overflow.                    */
long eatol(char *s)
{
    unsigned long n;
    int i;
    int sign = 1L;

    for (i = 0; isspace(s[i]); ++i) ;	/* skip white space */

    if (s[i] == '+' || s[i] == '-') {
	sign = (s[i++] == '+') ? 1L : -1L;
    }

    for (n = 0L; isdigit(s[i]) && n <= LONG_MAX/10L; ++i) {
	n = 10L * n + (s[i] - '0');
    }

    if ((sign ==  1L) && (n > LONG_MAX)   ||
	(sign == -1L) && (n > LONG_MIN) || isdigit(s[i]))
	sf_error("%s: eatol: overflow", __FILE__);

    return  sign * (long) n;
}


/* eatov - convert string s to unsigned long {0:ULONG_MAX}          *
 * Here, we check for overflow by seeing whether n decreases.       */
unsigned long eatov(char *s)
{
    unsigned long n;
    unsigned long n_old;
    int i;

    for (i = 0; isspace(s[i]); ++i) ;  /* skip white space */

    if (s[i] == '-')
	sf_error("%s: eatov: saw negative number", __FILE__);

    for (n_old = 0L, n = 0L; isdigit(s[i]); ++i) {
	n = 10L * n + (s[i] - '0');
	if (n < n_old)
	    sf_error("%s: eatov: overflow", __FILE__);
	n_old = n;
    }

    return n;
}

/* Dummy atof, atod routines until the ANSI police get here */
float eatof(char *s)
{
    return (float) atof(s);
}


double eatod(char *s)
{
    return atof(s);
}

#endif /* SUN_A */

#ifdef TEST
main(int argc, char **argv)
{
    char s[BUFSIZ];
    short nh;
    unsigned short nu;
    int ni;
    unsigned int np;
    long nl;
    unsigned long nv;

    initargs(argc, argv);


    /* Test code for eatoh */
    if (SHRT_MAX == LONG_MAX) {
	warn("Warning: eatoh not used on this machine.\n");
    } else {
	warn("\n");
    }
    strcpy(s, "0");
    nh = eatoh(s);
    warn("eatoh(%s) = %hd\n", s, nh);

    strcpy(s, "32767");
    nh = eatoh(s);
    warn("eatoh(%s) = %hd\n", s, nh);

    strcpy(s, "-32768");
    nh = eatoh(s);
    warn("eatoh(%s) = %hd\n", s, nh);


    /* Test code for eatou */
    if (USHRT_MAX == ULONG_MAX) {
	warn("Warning: eatou not used on this machine.\n");
    } else {
	warn("\n");
    }
    strcpy(s, "0");
    nu = eatou(s);
    warn("eatou(%s) = %hu\n", s, nu);

    strcpy(s, "65535");
    nu = eatou(s);
    warn("eatou(%s) = %hu\n", s, nu);


    /* Test code for eatoi */
    if (INT_MAX == LONG_MAX) {
	warn("Warning: eatoi not used on this machine.\n");
    } else {
	warn("\n");
    }
    strcpy(s, "0");
    ni = eatoi(s);
    warn("eatoi(%s) = %d\n", s, ni);

    strcpy(s, "2147483647");
    ni = eatoi(s);
    warn("eatoi(%s) = %d\n", s, ni);


    strcpy(s, "-2147483648");
    ni = eatoi(s);
    warn("eatoi(%s) = %d\n", s, ni);


    /* Test code for eatop */
    if (INT_MAX == LONG_MAX) {
	warn("Warning: eatop not used on this machine.\n");
    } else {
	warn("\n");
    }
    strcpy(s, "0");
    np = eatop(s);
    warn("eatop(%s) = %lu\n", s, np);

    strcpy(s, "4294967295");
    np = eatop(s);
    warn("eatop(%s) = %lu\n", s, np);


    /* Test code for eatol */
    warn("\n");
    strcpy(s, "0");
    nl = eatol(s);
    warn("eatol(%s) = %ld\n", s, nl);

    strcpy(s, "2147483647");
    nl = eatol(s);
    warn("eatol(%s) = %ld\n", s, nl);

    strcpy(s, "-2147483648");
    nl = eatol(s);
    warn("eatol(%s) = %ld\n", s, nl);


    /* Test code for eatov */
    strcpy(s, "0");
    nv = eatov(s);
    warn("eatov(%s) = %lu\n", s, nv);

    strcpy(s, "4294967295");
    nv = eatov(s);
    warn("eatov(%s) = %lu\n", s, nv);

    warn("Now we feed in 4294967296, expecting fatal error exit\n");
    strcpy(s, "4294967296");
    nv = eatov(s);
    warn("eatov(%s) = %lu\n", s, nv);

    return EXIT_SUCCESS;
}
#endif
