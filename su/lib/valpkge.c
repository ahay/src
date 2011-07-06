/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* VALPKGE: $Revision: 1.14 $ ; $Date: 1996/09/06 20:44:42 $	*/


/*********************** self documentation **********************/
/************************************************************************* 
VALPKGE - routines to handle variables of type Value

vtoi		cast Value variable as an int
vtol		cast Value variable as a long
vtof		cast Value variable as a float
vtod		cast Value variable as a double
atoval		convert ascii to Value
valtoabs	take absolute value of a Value variable
valcmp		compare Value variables
printfval	printf a Value variable
fprintfval	fprintf a Value variable
scanfval	scanf a Value variable
printftype	printf for the type of a segy header word

************************************************************************** 
Function Prototypes:
int vtoi(register cwp_String type, Value val);
long vtol(register cwp_String type, Value val);
float vtof(register cwp_String type, Value val);
double vtod(register cwp_String type, Value val);
void atoval(cwp_String type, cwp_String keyval, Value *valp);
Value valtoabs(cwp_String type, Value val);
int valcmp(register cwp_String type, Value val1, Value val2);
void printfval(register cwp_String type, Value val);
void fprintfval(FILE *stream, register cwp_String type, Value val);
void scanfval(register cwp_String type, Value *valp);

************************************************************************** 
Notes:
A Value is defined by the following in .../su/include/su.h:

typedef union { * storage for arbitrary type *
	char s[8];
	short h;
	unsigned short u;
	long l;
	unsigned long v;
	int i;
	unsigned int p;
	float f;
	double d;
	unsigned int U:16;
	unsigned int P:32;
} Value;

The use of the valpkge routines, as well as the hdrpkge routines,
permits the user to change the definition of the types of the 
various fields of the segy data type, without breaking codes that
look at part or all of these fields.

************************************************************************** 
Authors: CWP: Jack K. Cohen, Shuki Ronen
**************************************************************************/
/**************** end self doc ********************************/

#include <rsf.h>

#include "valpkge.h"
#include "atopkge.h"

#ifndef _sf_valpkge_h

/* TYPEDEFS */
typedef union { /* storage for arbitrary type */
	char s[8];
	short h;
	unsigned short u;
	long l;
	unsigned long v;
	int i;
	unsigned int p;
	float f;
	double d;
	unsigned int U:16;
	unsigned int P:32;
} Value;
/*^*/

typedef char *cwp_String;
/*^*/

#endif

int vtoi(register cwp_String type, Value val)
/*< cast Value variable as an int >*/
{
	switch(*type) {
		case 's': return (int) val.s[0];
		case 'h': return (int) val.h;
		case 'u': return (int) val.u;
		case 'l': return (int) val.l;
		case 'v': return (int) val.v;
		case 'i': return       val.i;
		case 'p': return (int) val.p;
		case 'f': return (int) val.f;
		case 'U': return (int) val.U;
		case 'P': return (int) val.P;
		default: sf_error("vtoi: unknown type %s", type);
			 return 0;	/* for lint */
	}
}

long vtol(register cwp_String type, Value val)
{
	switch(*type) {
		case 's': return (long) val.s[0];
		case 'h': return (long) val.h;
		case 'u': return (long) val.u;
		case 'l': return        val.l;
		case 'v': return (long) val.v;
		case 'i': return (long) val.i;
		case 'p': return (long) val.p;
		case 'f': return (long) val.f;
		case 'd': return (long) val.d;
		case 'U': return (long) val.U;
		case 'P': return (long) val.P;
		default: sf_error("vtol: unknown type %s", type);
			 return 0L;	/* for lint */
	}
}


float vtof(register cwp_String type, Value val)
/*< cast Value variable as a float >*/
{
	switch(*type) {
		case 's': return (float) val.s[0];
		case 'h': return (float) val.h;
		case 'u': return (float) val.u;
		case 'l': return (float) val.l;
		case 'v': return (float) val.v;
		case 'i': return (float) val.i;
		case 'p': return (float) val.p;
		case 'f': return         val.f;
		case 'd': return (float) val.d;
		case 'U': return (float) val.U;
		case 'P': return (float) val.P;
		default: sf_error("vtof: unknown type %s", type);
			 return 0.0;	/* for lint */
	}
}

double vtod(register cwp_String type, Value val)
{
	switch(*type) {
		case 's': return (double) val.s[0];
		case 'h': return (double) val.h;
		case 'u': return (double) val.u;
		case 'l': return (double) val.l;
		case 'v': return (double) val.v;
		case 'i': return (double) val.i;
		case 'p': return (double) val.p;
		case 'f': return (double) val.f;
		case 'd': return          val.d;
		case 'U': return (double) val.U;
		case 'P': return (double) val.P;
		default: sf_error("vtod: unknown type %s", type);
			 return 0.0;	/* for lint */
	}
}


int valcmp(register cwp_String type, Value val1, Value val2)
/*< compare Value variables >*/
{
	switch(*type) {
		case 's': return strcmp(val1.s, val2.s);
		case 'h':
			if      ( val1.h < val2.h ) return -1;
			else if ( val1.h > val2.h ) return  1;
			else                        return  0;
		case 'u':
			if      ( val1.u < val2.u ) return -1;
			else if ( val1.u > val2.u ) return  1;
			else                        return  0;
		case 'l':
			if      ( val1.l < val2.l ) return -1;
			else if ( val1.l > val2.l ) return  1;
			else                        return  0;
		case 'v':
			if      ( val1.v < val2.v ) return -1;
			else if ( val1.v > val2.v ) return  1;
			else                        return  0;
		case 'i':
			if      ( val1.i < val2.i ) return -1;
			else if ( val1.i > val2.i ) return  1;
			else                        return  0;
		case 'p':
			if      ( val1.p < val2.p ) return -1;
			else if ( val1.p > val2.p ) return  1;
			else                        return  0;
		case 'f':
			if      ( val1.f < val2.f ) return -1;
			else if ( val1.f > val2.f ) return  1;
			else                        return  0;
		case 'd':
			if      ( val1.d < val2.d ) return -1;
			else if ( val1.d > val2.d ) return  1;
			else                        return  0;
		case 'U':
			if      ( val1.U < val2.U ) return -1;
			else if ( val1.U > val2.U ) return  1;
			else                        return  0;
		case 'P':
			if      ( val1.P < val2.P ) return -1;
			else if ( val1.P > val2.P ) return  1;
			else                        return  0;
		default: sf_error("valcmp: unknown type %s", type);
			 return 0;	/* for lint */
	}
}


void printfval(register cwp_String type, Value val)
/*< printf a Value variable >*/
{
	switch(*type) {
	case 's':
		(void) printf("%s", val.s);
	break;
	case 'h':
		(void) printf("%d", val.h);
	break;
	case 'u':
		(void) printf("%d", val.u);
	break;
	case 'i':
		(void) printf("%d", val.i);
	break;
	case 'p':
		(void) printf("%d", val.p);
	break;
	case 'l':
		(void) printf("%ld", val.l);
	break;
	case 'v':
		(void) printf("%ld", val.v);
	break;
	case 'f':
		(void) printf("%f", val.f);
	break;
	case 'd':
		(void) printf("%f", val.d);
	break;
	case 'U':
		(void) printf("%d", val.U);
	break;
	case 'P':
		(void) printf("%d", val.P);
	break;
	default:
		sf_error("printfval: unknown type %s", type);
	}

	return;
}


void fprintfval(FILE *stream, register cwp_String type, Value val)
{
	switch(*type) {
	case 's':
		(void) fprintf(stream, "%s", val.s);
	break;
	case 'h':
		(void) fprintf(stream, "%d", val.h);
	break;
	case 'u':
		(void) fprintf(stream, "%d", val.u);
	break;
	case 'i':
		(void) fprintf(stream, "%d", val.i);
	break;
	case 'p':
		(void) fprintf(stream, "%d", val.p);
	break;
	case 'l':
		(void) fprintf(stream, "%ld", val.l);
	break;
	case 'v':
		(void) fprintf(stream, "%ld", val.v);
	break;
	case 'f':
		(void) fprintf(stream, "%f", val.f);
	break;
	case 'd':
		(void) fprintf(stream, "%f", val.d);
	break;
	case 'U':
		(void) fprintf(stream, "%d", val.U);
	break;
	case 'P':
		(void) fprintf(stream, "%d", val.P);
	break;
	default:
		sf_error("fprintfval: unknown type %s", type);
	}

	return;
}


void scanfval(register cwp_String type, Value *valp)
{
	switch(*type) {
	case 's':
		(void) scanf("%s", (char *) valp);
	break;
	case 'h':
	case 'u':
	case 'U':
		(void) scanf("%hd", (short *) valp);
	break;
	case 'i':
	case 'p':
	case 'P':
		(void) scanf("%d", (int *) valp);
	break;
	case 'l':
	case 'v':
		(void) scanf("%ld", (long int *)  valp);
	break;
	case 'f':
		(void) scanf("%f", (float *) valp);
	break;
	case 'd':
		(void) scanf("%lf", (double *) valp);
	break;
	default:
		sf_error("scanfval: unknown type %s", type);
	}

	return;
}




/* Convert ascii to Value according to type of keyword -- omitted bit types */
void atoval(
cwp_String type,	/* type of header keyword		*/
cwp_String keyval,	/* value of header keyword as ascii 	*/
Value *valp	/* pointer to converted value		*/
)
{
	switch(*type) {
	case 's':
		(void) strcpy(valp->s, keyval);
	break;
	case 'h':
		valp->h = eatoh(keyval); 
	break;
	case 'u':
		valp->u = eatou(keyval); 
	break;
	case 'i':
		valp->i = eatoi(keyval); 
	break;
	case 'p':
		valp->p = eatop(keyval); 
	break;
	case 'l':
		valp->l = eatol(keyval); 
	break;
	case 'v':
		valp->v = eatov(keyval); 
	break;
	case 'f':
		valp->f = eatof(keyval); 
	break;
	case 'd':
		valp->d = eatod(keyval); 
	break;
	default:
		sf_error("%s: %s: mysterious data type: %s",
					__FILE__, __LINE__, keyval);
	break;
	}

	return;
}

/* Get absolute value for type value variable */
Value valtoabs(cwp_String type, Value val)
{
        switch(*type) {
        case 'u':       /* do nothing if unsigned */
        case 'v':
        case 'p':
	case 'U':
	case 'P':
        break;
        case 'h':
                val.h = SF_ABS(val.h);
        break;
        case 'i':
                val.i = SF_ABS(val.i);
        break;
        case 'l':
                val.l = SF_ABS(val.l);
        break;
        case 'f':
                val.f = SF_ABS(val.f);
        break;
        case 'd':
                val.d = SF_ABS(val.d);
        break;
        default:
                sf_error("valtoabs: %d: mysterious type %s", __LINE__, type);
        }

        return val;
}

