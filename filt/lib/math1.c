#include <math.h>
#include <string.h>
#include <ctype.h>

#include "math1.h"
#include "stack.h"
#include "error.h"
#include "alloc.h"
#include "file.h"

typedef float (*func)(float);
static func functable[] = {
    cosf,
    sinf,
    tanf,
    acosf,
    asinf,
    atanf,
    expf,
    logf,
    sqrtf,
    fabsf
};

enum {GRP, NUM, INDX, FUN, UNARY, POW, MULDIV, PLUSMIN};

static sf_stack st1, st2;

static void check (void);

void sf_math_evaluate (int len, int nbuf, float** fbuf, float** fst)
{
    char *op;
    int *indx, i;
    float *num, f, *farr;
    func fun;

    sf_stack_set(st2,len);

    while (sf_full (st2)) {
	switch (sf_top (st2)) {
	    case NUM: 
		farr = *(++fst);
		num = (float*) sf_pop(st2);
		f = *num;
		for (i=0; i < nbuf; i++) { farr[i] = f; }
		break;
	    case INDX:
		farr = *(++fst);
		indx = (int*) sf_pop(st2);
		/* convert index to number */
		num = *(fbuf + (*indx));
		for (i=0; i < nbuf; i++) { farr[i] = num[i]; }
		break;
	    case FUN:
		indx = (int*) sf_pop(st2);
		fun = functable[*indx];
		farr = *fst;
		for (i=0; i < nbuf; i++) { farr[i] = fun(farr[i]); }
		break;
	    case UNARY:
		op = (char*) sf_pop(st2);
		farr = *fst;
		if ('-' == *op) {
		    for (i=0; i < nbuf; i++) { farr[i] = -farr[i]; }
		}
		break;
	    case POW:
		op = (char*) sf_pop(st2);
		num = *fst;
		farr = *(--fst);
		if ('^' == *op) {
		    for (i=0; i < nbuf; i++) { 
			farr[i] = powf(farr[i],num[i]);
		    }
		} else {
		    for (i=0; i < nbuf; i++) { 
			farr[i] = atan2f(farr[i],num[i]);
		    }
		}
		break;
	    case MULDIV:
		op = (char*) sf_pop(st2);
		num = *fst;
		farr = *(--fst);
		if ('*' == *op) {
		    for (i=0; i < nbuf; i++) { farr[i] *= num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { farr[i] /= num[i]; }
		}
		break;
	    case PLUSMIN:
		op = (char*) sf_pop(st2);
		num = *fst;
		farr = *(--fst);
		if ('+' == *op) {
		    for (i=0; i < nbuf; i++) { farr[i] += num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { farr[i] -= num[i]; }
		}
		break;
	    default:
		sf_error ("%s: syntax error in output",__FILE__);
		break;
	}
    }
}

int sf_math_parse (char* output, sf_file out)
{
    int i, j, keylen, *indx, type=-1, top, len, c, c2;
    char *key;
    float *num;
    bool hasleft;

    len = strlen(output);
    st1 = sf_stack_init (len);
    st2 = sf_stack_init (len);
    
    hasleft = false;
    for (i=0; i < len; i++) {
	c = output[i];
	
	if (isspace(c)) continue; /* skip white space */ 
	
	/* handle parentheses */

	if ('(' == c) {
	    hasleft = false;
	    sf_push(st2,&c,GRP);
	    continue;
	}

	if (')' == c) {
	    hasleft = true;
	    top = -1;
	    while (sf_full(st2)) {
		top = sf_top(st2);		
		if (GRP == top) {
		    sf_pop(st2);
		    break;
		}
		sf_push(st1,sf_pop(st2),top);
	    }
	    if (GRP != top) 
		sf_error("%s: unbalanced ')', position %d in output",
			 __FILE__,i);
	    continue;
	}
	
	if ('.' == c || isdigit(c)) { /* number */
	    hasleft = true;
	    for (j=i+1; j < len; j++) {
		c2 = output[j];
		if ('.' != c2 && !isdigit(c2)) break;
	    }
	    keylen = j-i;
	    key = sf_charalloc(keylen+1);
	    strncpy(key,output+i,keylen);
	    key[keylen]='\0';
 
	    num = sf_floatalloc(1);
	    if (0 == sscanf(key,"%f",num)) 
		sf_error("%s: ill-formated number: %s, position %d in output",
			 __FILE__,key,i);
	    free(key);

	    sf_push(st1,num,NUM); /* st1 is out */

	    i=j-1;
	    continue;
	}

	if (isalpha (c)) { /* identifier */
	    for (j=i+1; j < len; j++) {
		if (!isalnum((int) output[j])) break;
	    }

	    keylen = j-i;
	    key = sf_charalloc(keylen+1);
	    strncpy(key,output+i,keylen);
	    key[keylen]='\0';

	    indx = sf_intalloc(1);

	    if (sf_histint(out,key,indx)) {
		hasleft = true;
		sf_push(st1,indx,INDX); 
	    } else {
		hasleft = false;
		if (       0==strcmp(key,"cos"))  { *indx = 0;
		} else if (0==strcmp(key,"sin"))  { *indx = 1;
		} else if (0==strcmp(key,"tan"))  { *indx = 2;    
		} else if (0==strcmp(key,"acos")) { *indx = 3;
		} else if (0==strcmp(key,"asin")) { *indx = 4;
		} else if (0==strcmp(key,"atan")) { *indx = 5;
		} else if (0==strcmp(key,"exp"))  { *indx = 6;
		} else if (0==strcmp(key,"log"))  { *indx = 7;
		} else if (0==strcmp(key,"sqrt")) { *indx = 8;
		} else if (0==strcmp(key,"abs"))  { *indx = 9;
		} else {
		    sf_error("%s: unrecognized identifier: "
			     "%s, position %d in output",
			     __FILE__,key,i);
		}
		sf_push(st2,indx,FUN); 
	    }
	    free (key);
	    
	    i=j-1;
	    continue;
	}

	/* handle operators */
	switch (c) {
	    case '+': 
	    case '-':
		top = sf_top (st1);
		type = hasleft? PLUSMIN: UNARY;
		break;
	    case '*':
	    case '/':
		type = MULDIV;
		break;
	    case '^':
	    case '&':
		type = POW;
		break;
	    default:
		sf_error("%s: wrong character: %c, position %d in output",
			 __FILE__,c,i);
		break;
	}

	hasleft=false;

	while (sf_full(st2)) {
	    top = sf_top (st2);
	    if (GRP==top || top > type) break; /* compare precedence */
	    sf_push(st1,sf_pop(st2),top);
	}
	
	sf_push(st2,output+i,type);
    }
	
    /* push operators into output (reverse polish) */
    while (sf_full(st2)) {
	top = sf_top (st2);
	if (GRP == top)
	    sf_error("%s: unbalanced '(' in output",__FILE__);
	sf_push(st1,sf_pop(st2),top);
    }

    /* arrange in reverse order for computation */
    while (sf_full(st1)) {
	top = sf_top(st1);
	sf_push(st2,sf_pop(st1),top);
    }

    len = sf_stack_get(st2);
    check();
    sf_stack_close(st1);

    return len;
}

static void check (void)
{
    int top;

    while (sf_full (st2)) {
	switch (sf_top (st2)) {
	    case NUM: 
	    case INDX:
		sf_push(st1,sf_pop(st2),NUM);
		break;
	    case FUN:
	    case UNARY:
		sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("%s [%d]: syntax error in output",
			      __FILE__,__LINE__);
		break;
	    case POW:
	    case MULDIV:
	    case PLUSMIN:
		sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("%s[%d]: syntax error in output",
			      __FILE__,__LINE__);
		sf_pop(st1);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("%s[%d]: syntax error in output",
			      __FILE__,__LINE__);
		break;
	    default:
		sf_error ("%s[%d]: syntax error in output",
			  __FILE__,__LINE__);
		break;
	}
    }
    top = sf_top (st1);
    if (0 != sf_stack_get(st1) || (NUM != top && INDX != top)) 
	sf_error ("%s[%d]: syntax error in output",
		  __FILE__,__LINE__);
}

/* 	$Id$	 */
