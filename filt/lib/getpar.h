#ifndef _sf_getpar_h
#define _sf_getpar_h

#include <stdio.h>

#include "c99.h"

void sf_init(int argc,char *argv[]);
void sf_parclose (void);
char* sf_getprog (void);
char* sf_getuser (void);
char* sf_gethost (void);
bool sf_getint (const char* key,/*@out@*/ int* par);
bool sf_getints (const char* key,/*@out@*/ int* par,size_t n);
bool sf_getfloat (const char* key,/*@out@*/ float* par);
bool sf_getfloats (const char* key,/*@out@*/ float* par,size_t n);
char* sf_getstring (const char* key);
bool sf_getstrings (const char* key,/*@out@*/ char** par,size_t n);
bool sf_getbool (const char* key,/*@out@*/ bool* par);
bool sf_getbools (const char* key,/*@out@*/ bool* par,size_t n);

#endif

/* 	$Id: getpar.h,v 1.5 2004/03/22 05:43:24 fomels Exp $	 */
