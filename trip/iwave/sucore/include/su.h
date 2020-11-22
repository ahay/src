/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/* su.h - include file for SU programs
 *
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/su.h,v $
 * $Revision: 1.33 $ ; $Date: 1997/10/15 15:18:21 $
 */

#ifndef SU_H
#define SU_H

#include "par.h"

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


/* DEFINES */
#define CHECK_NT(label,nt) \
	if(nt > SU_NFLTS) suerr("%s=%d must not exceed %d",label,nt,SU_NFLTS)
#define NALLOC	(524288)
#define NFALLOC	(NALLOC/FSIZE)
#define NIALLOC	(NALLOC/ISIZE)
#define NDALLOC	(NALLOC/DSIZE)
#define LOWBYTE(w) ((w) & 0xFF)
#define HIGHBYTE(w) LOWBYTE((w) >>8)
#define LOWWORD(w) ((w) & 0xFFFF)
#define HIGHWORD(w) LOWWORD((w) >>16)
#define ISNEGCHAR(c) ((c) & 0x80)
#define SIGNEXTEND(c) (~0xFF | (int) (c))

/*	READ_OK  - read  permission for access(2)
 *	WRITE_OK - write permission for access(2)
 *	EXEC_OK  - exec  permission for access(2)
 *	FILE_OK  - file  existence  for access(2)
 *	Note: these are changed from the usual defines in file.h
 *	      because this include exists on some machines and
 *	      not others, often overlaps fcntl.h, etc.  Lint is
 *            happier with a fresh start.
 *	Note: Post-ANSI sometimes R_OK in unistd.h (this isn't
 *	      an ANSI file).
 */
#define		READ_OK		4
#define		WRITE_OK	2
#define		EXEC_OK		1
#define		FILE_OK		0

/* For plotting by keyword */
#define	 IS_DEPTH(str)	((  STREQ(str,"gelev")	|| \
			   STREQ(str,"selev")	|| \
			   STREQ(str,"sdepth")	|| \
			   STREQ(str,"gdel")	|| \
			   STREQ(str,"sdel")	|| \
			   STREQ(str,"swdep")	|| \
			   STREQ(str,"gwdep")	      )?cwp_true:cwp_false)

#define	 IS_COORD(str)	(( STREQ(str,"sx")   || \
			   STREQ(str,"sy")   || \
			   STREQ(str,"gx")   || \
			   STREQ(str,"gy")	   )?cwp_true:cwp_false)

/* FUNCTION PROTOTYPES */
#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/* valpkge */
int vtoi(cwp_String type, Value val);
long vtol(cwp_String type, Value val);
float vtof(cwp_String type, Value val);
double vtod(cwp_String type, Value val);
int valcmp(cwp_String type, Value val1, Value val2);
void printfval(cwp_String type, Value val);
void fprintfval(FILE *stream, cwp_String type, Value val);
void scanfval(cwp_String type, Value *valp);
void atoval(cwp_String type, cwp_String keyval, Value *valp);
void getparval(cwp_String name, cwp_String type, int n, Value *valp);
Value valtoabs(cwp_String type, Value val);

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif
