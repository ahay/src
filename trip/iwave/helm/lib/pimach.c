/* pimach.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

double pimach_(float *dum)
{
    /* System generated locals */
    float ret_val;

/* ***BEGIN PROLOGUE  PIMACH */

/*     This subprogram supplies the value of the constant PI correct to */
/*     machine precision where */

/*     PI=3.1415926535897932384626433832795028841971693993751058209749446 */
/* ***ROUTINES CALLED  (NONE) */


/* ***END PROLOGUE  PIMACH */

/* ***FIRST EXECUTABLE STATEMENT  PIMACH */
    ret_val = 3.14159265358979f;
    return ret_val;
} /* pimach_ */

