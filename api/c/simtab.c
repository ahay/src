/* Simbol Table for parameters. Implemented as a hash table. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <ctype.h>

#include <stdio.h>
/*^*/

#include "simtab.h"
#include "alloc.h"
#include "error.h"

#include "_bool.h"
/*^*/

#ifndef _sf_simtab_h

#define SF_EOL '\014'
#define SF_EOT '\004'
/*^*/

typedef struct sf_SimTab *sf_simtab; /* Simbol Table structure */
/*^*/

#endif

#ifdef _POSIX_ARG_MAX
#define LINELEN _POSIX_ARG_MAX
#else
#ifdef ARG_MAX
#define LINELEN ARG_MAX
#else
#define LINELEN 1024
#endif
#endif

struct entry {
    char *key, *val;
    struct entry* next;
};

struct sf_SimTab {
    struct entry **pars;
    int size;
};


static unsigned int hash (const char *key, int size)
/* Taken from Kernigan and Pike, "The practice of programming" */
{
    unsigned int h;
    unsigned char *p;

    h=0;
    for (p = (unsigned char*) key; *p != '\0'; p++) {
	h = 31 * h + (int) *p;
    }
    return (h % size);
}

static char *
strsep1(char **stringp, const char *delim)
/* portable strsep */
{
    char *start = *stringp;
    char *cp;
    char ch;

    if (start == NULL) return NULL;

    for (cp = start; (ch = *cp); cp++) {
        if (strchr(delim, ch)) {
            *cp++ = 0;
            *stringp = cp;
            return start;
        }
    }
    *stringp = NULL;
    return start;
}

sf_simtab sf_simtab_init(int size)
/*< Create simbol table. >*/
{
    sf_simtab table;
    int i;


    table = (sf_simtab) sf_alloc(1,sizeof(*table));
    table->pars = (struct entry **) sf_alloc (size,sizeof(struct entry*));
    table->size = size;

    for (i=0; i < size; i++) {
	table->pars[i] = NULL;
    }

    return table;
}

void sf_simtab_close(sf_simtab table)
/*< Free allocated memory >*/
{
    int i;
    struct entry *e, *next;

    if (NULL==table) return;

    for (i=0; i < table->size; i++) {
	for (e = table->pars[i]; e != NULL; e = next) {	    
	    free (e->key);
	    free (e->val);
	    next = e->next;
	    free (e);
	}	
    }
    free(table->pars);
    free(table);
}

void sf_simtab_enter(sf_simtab table, const char *key, const char* val)
/*< Add an entry key=val to the table >*/
{
    int h;
    struct entry *e;
    size_t len;

    h = hash(key,table->size);

    /* check if exists */
    
    for (e = table->pars[h]; e != NULL; e = e->next) {
	if (strcmp(key,e->key) == 0) {
	    if (strcmp(val,e->val) != 0) {
		len = strlen(val)+1;
		e->val = (char*) sf_realloc(e->val,len,1);
		memcpy(e->val,val,len);
	    }
	    return;
	}
    }

    /* create new */
    
    e = (struct entry*) sf_alloc (1,sizeof(*e));

    len = strlen(key)+1;
    e->key = sf_charalloc(len); 
    memcpy(e->key,key,len);
    
    len = strlen(val)+1;
    e->val = sf_charalloc(len); 
    memcpy(e->val,val,len);

    e->next = table->pars[h];
    table->pars[h] = e;
}

char *sf_simtab_get(sf_simtab table, const char *key) 
/*< extract a value from the table >*/
{
    int h;
    struct entry *e;

    h = hash(key,table->size);

    for (e = table->pars[h]; e != NULL; e = e->next) {
	if (0 == strcmp(key,e->key)) return e->val;
    }
    return NULL;
}

bool sf_simtab_getint (sf_simtab table, const char* key,/*@out@*/ int* par)
/*< extract an int parameter from the table >*/
{
    char* val;
    long int i;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;
 
    errno = 0;
    i = strtol(val,NULL,10);
    if (ERANGE == errno || i < INT_MIN || i > INT_MAX) return false;

    *par = (int) i;
    return true;
}

bool sf_simtab_getlargeint (sf_simtab table, const char* key,/*@out@*/ off_t* par)
/*< extract a sf_largeint parameter from the table >*/
{
    char* val;
    off_t i;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;
    
    errno = 0;
    i = strtoll(val,NULL,10);
    if (ERANGE == errno || i < LONG_MIN || i > LONG_MAX) return false;

    *par = (off_t) i;
    return true;
}

bool sf_simtab_getfloat (sf_simtab table, const char* key,/*@out@*/ float* par)
/*< extract a float parameter from the table >*/
{
    char* val;
    double f;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;

    errno = 0;
    f = strtod(val,NULL);
    if (ERANGE == errno || f < -FLT_MAX || f > FLT_MAX) return false;
    
    *par = (float) f;
    return true;
}

bool sf_simtab_getdouble (sf_simtab table, const char* key,/*@out@*/ double* par)
/*< extract a double parameter from the table >*/
{
    char* val;
    double f;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;

    errno = 0;
    f = strtod(val,NULL);
    if (ERANGE == errno) return false;
    
    *par = f;
    return true;
}

bool sf_simtab_getfloats (sf_simtab table, const char* key,
			  /*@out@*/ float* par,size_t n)
/*< extract a float array parameter from the table >*/
{
    size_t i, clen;
    long num;
    char *val, *fval, *cnum, *fvali;
    double fi=0.;
    bool result;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;
  
    cnum = sf_charalloc(strlen(val));
    result = true;
    for (i = 0; i < n; i++) {
	fval = (0==i)? strtok(val,","):strtok(NULL,",");
	if (NULL == fval) {
	    if (0==i) return false;
	} else {
	    fvali = strpbrk(fval,"x*");
	    if (NULL != fvali) {
		clen = fvali-fval;
		strncpy(cnum,fval,clen);
		cnum[clen]='\0';

		num = strtol(cnum,NULL,10);
		if (ERANGE == errno) 
		    sf_error ("%s: Wrong counter in %s='%s':",__FILE__,
			      key,fval);
		fvali++;

		fi = strtod(fvali,NULL);
		if (ERANGE == errno || fi < -FLT_MAX || fi > FLT_MAX ) 
		    sf_error ("%s: %s='%s' is out of range:",__FILE__,
			      key,fvali);
		
		for (; i < n && i < (size_t) num; i++) {
		    par[i] = (float) fi;
		}
		if (i == n) break;
	    } else {
		if (NULL==fval) {
		    if (0==i) {
			result = false;
			break;
		    }
		} else {
		    fi = strtod(fval,NULL);
		    if (ERANGE == errno || fi < -FLT_MAX || fi > FLT_MAX ) 
			sf_error("%s: %s='%s' is out of range:",
				 __FILE__,key,fval);
		}
	    }
	}
	par[i] = (float) fi;
    }
    free(cnum);

    return result;
}

char* sf_simtab_getstring (sf_simtab table, const char* key) 
/*< extract a string parameter from the table >*/
{
    char *val, *string, *qopen, *qclose;
    int iopen;
    size_t iclose;

    val = sf_simtab_get(table,key);
    if (NULL == val) return NULL;

    qopen = strchr(val,'\"');   
    iopen = (NULL == qopen)? 0: (qopen-val)+1; 	
    qclose = strrchr(val,'\"'); 
    iclose = (NULL == qclose)? strlen(val): (size_t) (qclose-val);
    iclose -= iopen;
    string = sf_charalloc (iclose+1);
    memcpy(string,val+iopen,iclose);
    string[iclose]='\0';

    return string;
}

bool sf_simtab_getbool (sf_simtab table, const char* key,/*@out@*/ bool *par)
/*< extract a bool parameter from the table >*/
{
    char* val;

    val = sf_simtab_getstring(table,key);
    if (NULL == val) return false;

    *par = (bool) (('y' == val[0]) || 
		   ('Y' == val[0]) || 
		   ('1' == val[0]));

    free(val);
    return true;
}

bool sf_simtab_getbools (sf_simtab table, const char* key,
			 /*@out@*/bool *par,size_t n)
/*< extract a bool array parameter from the table >*/
{
    bool test=false;
    size_t i, clen;
    long num;
    char *val, *cnum, *fval, *fvali;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;

    cnum = sf_charalloc(strlen(val));
    for (i = 0; i < n; i++) {
	fval = (0==i)? strtok(val,","):strtok(NULL,",");
	fvali = strpbrk(fval,"x*");
	if (NULL != fvali) {
	    clen = fvali-fval;
	    strncpy(cnum,fval,clen);
	    cnum[clen]='\0';

	    num = strtol(cnum,NULL,10);
	    if (ERANGE == errno) 
		sf_error ("%s: Wrong counter in %s='%s':",__FILE__,key,fval);
	    fvali++;
	    test = (bool) (('y' == fvali[0]) || 
			   ('Y' == fvali[0]) || 
			   ('1' == fvali[0]));
	    for (; i < n && i < (size_t) num; i++) {
		par[i] = test;
	    }
	    if (i == n) return true;
	} else {
	    if (NULL==fval) {
		if (0==i) return false;
	    } else {
		test = (bool) (('y' == fval[0]) || 
			       ('Y' == fval[0]) || 
			       ('1' == fval[0]));
	    }
	}
	par[i] = test;
    }
    free(cnum);

    return true;
}

bool sf_simtab_getints (sf_simtab table, const char* key,
			/*@out@*/ int *par,size_t n)
/*< extract an int array parameter from the table >*/
{    
    size_t i, clen;
    long num, j=0;
    char *val, *cnum, *fval, *fvali;

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;
    
    cnum = sf_charalloc(strlen(val));
    for (i = 0; i < n; i++) {
	fval = (0==i)? strtok(val,","):strtok(NULL,",");
	if (NULL==fval) {
	    if (0==i) return false;
	} else {
	    fvali = strpbrk(fval,"x*");
	    if (NULL != fvali) {
		clen = fvali-fval;
		strncpy(cnum,fval,clen);
		cnum[clen]='\0';

		num = strtol(cnum,NULL,10);
		if (ERANGE == errno) 
		    sf_error ("%s: Wrong counter in %s='%s':",
			      __FILE__,key,fval);
		fvali++;
		j = strtol(fvali,NULL,10);
		if (ERANGE == errno || j < INT_MIN || j > INT_MAX) 
		    sf_error ("%s: Wrong value in %s='%s':",__FILE__,key,fval);
		for (; i < n && i < (size_t) num; i++) {
		    par[i] = (int) j;
		}
		if (i == n) {
		    free(cnum);
		    return true;
		}
	    } else {
		j = strtol(fval,NULL,10);
		if (ERANGE == errno || j < INT_MIN || j > INT_MAX) 
		    sf_error ("%s: Wrong value in %s='%s':",__FILE__,key,fval);
	    }
	}
	par[i] = (int) j;
    }
    free(cnum);

    return true;
}

bool sf_simtab_getstrings (sf_simtab table, const char* key,
			   /*@out@*/ char **par,size_t n)
/*< extract a string array parameter from the table >*/
{    
    int i, iopen;
    size_t iclose;
    char *val, *string, *qopen, *qclose; 

    val = sf_simtab_get(table,key);
    if (NULL == val) return false;

    qopen = strchr(val,'\"');   
    iopen = (NULL == qopen)? 0: (qopen-val)+1; 	
    qclose = strrchr(val,'\"'); 
    iclose = (NULL == qclose)? strlen(val): (size_t) (qclose-val);
    iclose -= iopen;
    string = sf_charalloc (iclose+1);
    memcpy(string,val+iopen,iclose);
    string[iclose]='\0';
   
    par[0] = strsep1(&string,":");
    for (i = 1; i < n; i++) {
	par[i] = strsep1(&string,":");
    }

    return true;
}

void sf_simtab_put (sf_simtab table, const char *keyval) 
/*< put a key=val string to the table >*/
{
    char *eq, *key;
    size_t keylen;

    eq = strchr(keyval,'=');
    if (NULL == eq) return;
    eq++;
    
    keylen = (size_t) (eq-keyval);
    key = sf_charalloc(keylen);
    memcpy(key,keyval,keylen);
    key[keylen-1]='\0';

    sf_simtab_enter(table,key,eq);
    free(key);
}

void sf_simtab_string (sf_simtab table, char* string)
/*< extract parameters from a string >*/
{
    char *cw, *cl, word[LINELEN];
    int c;
    enum {START, INAWORD, STRING} state;

    cw = word;
    state = START;
    for (cl = string; '\0' != (c=*cl); cl++) {
	switch (state) {
	    case START:
		if (!isspace(c)) {
		    *cw++ = c;
		    state = ('"' == (char) c)? STRING:INAWORD;
		}
		break;
	    case INAWORD:
		if (isspace(c)) {
		    *cw = '\0';
		    sf_simtab_put (table,word);
		    cw = word;
		    state = START;
		} else {
		    *cw++ = c;
		    if ('"' == c) state = STRING;
		}
		break;
	    case STRING:
		*cw++ = c;
		if ('"' == c) state = INAWORD;
		break;
	}
    }
     
    if (STRING == state) {
	sf_error("%s: unterminated string",__FILE__);
    } else if (INAWORD == state) {
	*cw = '\0';
	sf_simtab_put (table,word);
    }
}

void sf_simtab_input (sf_simtab table, FILE* fp, FILE* out) 
/*< extract parameters from a file >*/
{
    char line[LINELEN];

    while (NULL != fgets(line,4,fp)) {
	/* code for the header end */
	if (SF_EOL==line[0] && SF_EOL==line[1] && SF_EOT==line[2]) break;
	if ('\n'  !=line[0] && '\n'  !=line[1] && '\n'  !=line[2] &&
	    NULL == fgets(line+3,LINELEN-3,fp)) break;

	if (NULL != out) fputs(line,out);
	sf_simtab_string(table,line);
    }
    if (NULL != out) fflush(out);
}

void sf_simtab_output (sf_simtab table, FILE* fp) 
/*< output parameters to a file >*/
{
    int i;
    struct entry *e;

    for (i=0; i < table->size; i++) {
	for (e = table->pars[i]; e != NULL; e = e->next) {	  
	    if (0 >= fprintf(fp,"\t%s=%s\n",e->key,e->val))
		sf_error ("%s: Trouble printing table:",__FILE__);
	}
    }
}

void sf_simtab_expand (sf_simtab table, const sf_simtab expand) 
/*< expand table from another table >*/
{
    int i;
    struct entry *e;

    for (i=0; i < expand->size; i++) {
	for (e = expand->pars[i]; e != NULL; e = e->next) {
	    sf_simtab_enter(table, e->key, e->val);
	}
    }
}

/* 	$Id$	 */
