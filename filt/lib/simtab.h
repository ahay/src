#ifndef _sf_simtab_h
#define _sf_simtab_h

#include <stdio.h>

#include "c99.h"

/* Simbol Table structure */
typedef struct sf_SimTab *sf_simtab;

sf_simtab sf_simtab_init(int n);
void sf_simtab_close(sf_simtab table);
char *sf_simtab_get(sf_simtab table, const char *key);
bool sf_simtab_getint (sf_simtab table, const char* key,/*@out@*/ int* par);
bool sf_simtab_getints (sf_simtab table, const char* key,
			/*@out@*/ int* par,size_t n);
bool sf_simtab_getfloat (sf_simtab table, const char* key,/*@out@*/ float* par);
bool sf_simtab_getfloats (sf_simtab table, const char* key,
			  /*@out@*/ float* par,size_t n);
bool sf_simtab_getbool (sf_simtab table, const char* key,/*@out@*/ bool* par);
bool sf_simtab_getbools (sf_simtab table, const char* key,
			 /*@out@*/ bool* par,size_t n);
char* sf_simtab_getstring (sf_simtab table, const char* key);
bool sf_simtab_getstrings (sf_simtab table, const char* key,
			   /*@out@*/ char** par,size_t n);
void sf_simtab_enter(sf_simtab table, const char *key, const char* val);
void sf_simtab_put(sf_simtab table, const char *keyval);
void sf_simtab_input (sf_simtab table, FILE* fp);
void sf_simtab_output (sf_simtab table, FILE* fp);

#endif

