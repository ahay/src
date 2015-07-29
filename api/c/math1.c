/* Evaluating mathematical expressions. */
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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "math1.h"
#include "c99.h"
#include "stack.h"
#include "error.h"
#include "alloc.h"
#include "komplex.h"
#include "_defs.h"

#include "kiss_fft.h"
#include "file.h"
/*^*/

static float mysign(float f)
{
    f = (f >= 0  ?  1.0f  :  -1.0f );
    return f;
}

typedef float (*func)(float);
static func functable[] = {
    cosf,
    sinf,
    tanf,
    acosf,
    asinf,
    atanf,
    coshf,
    sinhf,
    tanhf,
    acoshf,
    asinhf,
    atanhf,
    expf,
    logf,
    sqrtf,
    fabsf,
    erff,
    erfcf,
    mysign
};

static sf_complex myabs(sf_complex c)
{
    c = sf_cmplx(cabsf(c),0.);
    return c;
}

static sf_complex myconj(sf_complex c)
{
    c = conjf(c);
    return c;
}

static sf_complex myarg(sf_complex c)
{
    c = sf_cmplx(cargf(c),0.);
    return c;
}

static sf_complex myreal(sf_complex c)
{
    c = sf_cmplx(crealf(c),0.);
    return c;
}

static sf_complex myimag(sf_complex c)
{
    c = sf_cmplx(0.,cimagf(c));
    return c;
}

typedef sf_complex (*cfunc)(sf_complex);
static cfunc cfunctable[] = {
    ccosf,
    csinf,
    ctanf,
    cacosf,
    casinf,
    catanf,
    ccoshf,
    csinhf,
    ctanhf,
    cacoshf,
    casinhf,
    catanhf,
    cexpf,
    clogf,
    csqrtf,
    myabs,
    myarg,
    myconj,
    myreal,
    myimag
/*    cerf - in the future */
};


static int myisign(int i)
{
    return SF_SIG(i);
}

static int myiabs(int i)
{
    return SF_ABS(i);
}

typedef int (*ifunc)(int);
static ifunc ifunctable[] = {
    myisign,
    myiabs
};
    
enum {GRP, NUM, INDX, FUN, POW, MOD, MULDIV, UNARY, PLUSMIN};

static sf_stack st1, st2;

static void check (void);


void sf_math_evaluate (int     len  /* stack length */, 
		       int     nbuf /* buffer length */, 
		       float** fbuf /* number buffers */, 
		       float** fst  /* stack */)
/*< Evaluate a mathematical expression from stack (float numbers) >*/
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

void sf_int_math_evaluate (int   len  /* stack length */, 
			   int   nbuf /* buffer length */, 
			   int** ibuf /* number buffers */, 
			   int** ist  /* stack */)
/*< Evaluate a mathematical expression from stack (integer numbers) >*/
{
    char *op;
    int *indx, i, *num, j, k, *iarr;
    ifunc fun;

    sf_stack_set(st2,len);

    while (sf_full (st2)) {
	switch (sf_top (st2)) {
	    case NUM: 
		iarr = *(++ist);
		num = (int*) sf_pop(st2);
		j = *num;
		for (i=0; i < nbuf; i++) { iarr[i] = j; }
		break;
	    case INDX:
		iarr = *(++ist);
		indx = (int*) sf_pop(st2);
		/* convert index to number */
		num = *(ibuf + (*indx));
		for (i=0; i < nbuf; i++) { iarr[i] = num[i]; }
		break;
	    case FUN:
		indx = (int*) sf_pop(st2);
		fun = ifunctable[*indx];
		iarr = *ist;
		for (i=0; i < nbuf; i++) { iarr[i] = fun(iarr[i]); }
		break;
	    case UNARY:
		op = (char*) sf_pop(st2);
		iarr = *ist;
		if ('-' == *op) {
		    for (i=0; i < nbuf; i++) { iarr[i] = -iarr[i]; }
		}
		break;
	    case POW:
		op = (char*) sf_pop(st2);
		num = *ist;
		iarr = *(--ist);
		if ('^' == *op) {
		    for (i=0; i < nbuf; i++) {
			j=1;
			for (k=0; k < num[i]; k++) {
			    j *= iarr[i];
			}
			iarr[i] = j;
		    }
		} else {
		    sf_error("unsupported \"%s\" operation",op);
		}
		break;
	    case MOD:
		op = (char*) sf_pop(st2);
		num = *ist;
		iarr = *(--ist);
		for (i=0; i < nbuf; i++) {
		    iarr[i] = iarr[i]%num[i];
		}
		break;
	    case MULDIV:
		op = (char*) sf_pop(st2);
		num = *ist;
		iarr = *(--ist);
		if ('*' == *op) {
		    for (i=0; i < nbuf; i++) { iarr[i] *= num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { iarr[i] /= num[i]; }
		}
		break;
	    case PLUSMIN:
		op = (char*) sf_pop(st2);
		num = *ist;
		iarr = *(--ist);
		if ('+' == *op) {
		    for (i=0; i < nbuf; i++) { iarr[i] += num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { iarr[i] -= num[i]; }
		}
		break;
	    default:
		sf_error ("%s: syntax error in output",__FILE__);
		break;
	}
    }
}


void sf_complex_math_evaluate (int          len  /* stack length */, 
			       int          nbuf /* buffer length */, 
			       sf_complex** cbuf /* number buffers */, 
			       sf_complex** cst  /* stack */)
/*< Evaluate a mathematical expression from stack (complex numbers) >*/
{
    char *op;
    int *indx, i;
    float *fnum;
    sf_complex *num, c, *carr;
    cfunc fun;
    
    sf_stack_set(st2,len);

    while (sf_full (st2)) {
	switch (sf_top (st2)) {
	    case NUM: 
		carr = *(++cst);
		fnum = (float*) sf_pop(st2);
		c = sf_cmplx(*fnum,0.);
		for (i=0; i < nbuf; i++) { carr[i] = c; }
		break;
	    case INDX:
		carr = *(++cst);
		indx = (int*) sf_pop(st2);
		/* convert index to number */
		num = *(cbuf + (*indx));
		for (i=0; i < nbuf; i++) { carr[i] = num[i]; }
		break;
	    case FUN:
		indx = (int*) sf_pop(st2);
		fun = cfunctable[*indx];
		carr = *cst;
		for (i=0; i < nbuf; i++) { carr[i] = fun(carr[i]); }
		break;
	    case UNARY:
		op = (char*) sf_pop(st2);
		carr = *cst;
		if ('-' == *op) {
		    for (i=0; i < nbuf; i++) {
#ifdef SF_HAS_COMPLEX_H 
			carr[i] = -carr[i];
#else
			carr[i] = sf_cneg(carr[i]);
#endif
		    }
		}
		break;
	    case POW:
		op = (char*) sf_pop(st2);
		num = *cst;
		carr = *(--cst);
		if ('^' == *op) {
		    for (i=0; i < nbuf; i++) { 
			carr[i] = cpowf(carr[i],num[i]);
		    }
		} 
		break;
	    case MULDIV:
		op = (char*) sf_pop(st2);
		num = *cst;
		carr = *(--cst);
		if ('*' == *op) {
		    for (i=0; i < nbuf; i++) { 
#ifdef SF_HAS_COMPLEX_H
			carr[i] *= num[i];
#else
			carr[i] = sf_cmul(carr[i],num[i]);
#endif 
		    }
		} else {
		    for (i=0; i < nbuf; i++) { 
#ifdef SF_HAS_COMPLEX_H
			carr[i] /= num[i]; 
#else
			carr[i] = sf_cdiv(carr[i],num[i]);
#endif
		    }
		}
		break;
	    case PLUSMIN:
		op = (char*) sf_pop(st2);
		num = *cst;
		carr = *(--cst);
		if ('+' == *op) {
		    for (i=0; i < nbuf; i++) {
#ifdef SF_HAS_COMPLEX_H
			carr[i] += num[i]; 
#else
			carr[i] = sf_cadd(carr[i],num[i]);
#endif
		    }
		} else {
		    for (i=0; i < nbuf; i++) { 
#ifdef SF_HAS_COMPLEX_H
			carr[i] -= num[i]; 
#else
			carr[i] = sf_csub(carr[i],num[i]);
#endif
		    }
		}
		break;
	    default:
		sf_error ("%s: syntax error in output",__FILE__);
		break;
	}
    }
}

size_t sf_math_parse (char*       output /* expression */, 
		      sf_file     out    /* parameter file */,
		      sf_datatype datatype)
/*< Parse a mathematical expression, returns stack length >*/ 
{
    int *indx, type=-1, top, c, c2, *inum;
    size_t i, j, keylen, len;
    char *key;
    float *num;
    bool hasleft;

    len = strlen(output);
    st1 = sf_stack_init (len);
    st2 = sf_stack_init (len);
    
    hasleft = false;
    for (i=0; i < len; i++) {
	c = (int) output[i];
	
	if (isspace(c)) continue; /* skip white space */ 
	
	/* handle parentheses */

	if ('(' == (char) c) {
	    hasleft = false;
	    sf_push(st2,&c,GRP);
	    continue;
	}

	if (')' == (char) c) {
	    hasleft = true;
	    top = -1;
	    while (sf_full(st2)) {
		top = sf_top(st2);		
		if (GRP == top) {
		    (void) sf_pop(st2);
		    break;
		}
		sf_push(st1,sf_pop(st2),top);
	    }
	    if (GRP != top) 
		sf_error("%s: unbalanced ')', position %d in output",
			 __FILE__,i);
	    continue;
	}
	
	if (SF_INT == datatype) {
	    if (isdigit(c)) { /* number */
		hasleft = true;
		c2 = ' ';
		for (j=i+1; j < len; j++) {
		    c2 = output[j];
		    if (!isdigit(c2)) break;
		}

		keylen = j-i;
		key = sf_charalloc(keylen+1);
		strncpy(key,output+i,keylen);
		key[keylen]='\0';
		
		inum = sf_intalloc(1);
		if (0 == sscanf(key,"%d",inum)) 
		    sf_error("%s: ill-formated number: %s, position %d in output",
			     __FILE__,key,i);
		free(key);
		
		sf_push(st1,inum,NUM); /* st1 is out */
		
		i=j-1;
		continue;
	    }
	} else {
	    if ('.' == (char) c || isdigit(c)) { /* number */
		hasleft = true;
		c2 = ' ';
		for (j=i+1; j < len; j++) {
		    c2 = output[j];
		    if ('.' != c2 && !isdigit(c2) && 
			'e' != c2) break;
		}
		if ('e' == output[j-1] && ('-' == c2 || '+' == c2)) {
		    for (j++; j < len; j++) {
			c2 = output[j];
			if (!isdigit(c2)) break;
		    }
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
		if (SF_INT==datatype) {
		    if        (0==strcmp(key,"sign")) { *indx = 0;
		    } else if (0==strcmp(key,"abs"))  { *indx = 1;
		    } else {
			sf_error("%s: unrecognized identifier: "
				 "%s, position %d in output",
				 __FILE__,key,i);
		    }
		} else if (0==strcmp(key,"cos"))   { *indx =  0;
		} else if (0==strcmp(key,"sin"))   { *indx =  1;
		} else if (0==strcmp(key,"tan"))   { *indx =  2;    
		} else if (0==strcmp(key,"acos"))  { *indx =  3;
		} else if (0==strcmp(key,"asin"))  { *indx =  4;
		} else if (0==strcmp(key,"atan"))  { *indx =  5;
		} else if (0==strcmp(key,"cosh"))  { *indx =  6;
		} else if (0==strcmp(key,"sinh"))  { *indx =  7;
		} else if (0==strcmp(key,"tanh"))  { *indx =  8;    
		} else if (0==strcmp(key,"acosh")) { *indx =  9;
		} else if (0==strcmp(key,"asinh")) { *indx = 10;
		} else if (0==strcmp(key,"atanh")) { *indx = 11;
		} else if (0==strcmp(key,"exp"))   { *indx = 12;
		} else if (0==strcmp(key,"log"))   { *indx = 13;
		} else if (0==strcmp(key,"sqrt"))  { *indx = 14;
		} else if (0==strcmp(key,"abs"))   { *indx = 15;
		} else if (0==strcmp(key,"erf")  && SF_FLOAT==datatype)   { *indx = 16;
		} else if (0==strcmp(key,"erfc") && SF_FLOAT==datatype)   { *indx = 17;
		} else if (0==strcmp(key,"sign") && SF_FLOAT==datatype)   { *indx = 18;   
		} else if (0==strcmp(key,"arg") && SF_COMPLEX==datatype) { *indx = 16;
		} else if (0==strcmp(key,"conj") && SF_COMPLEX==datatype) { *indx = 17;
		} else if (0==strcmp(key,"real") && SF_COMPLEX==datatype) { *indx = 18;
		} else if (0==strcmp(key,"imag") && SF_COMPLEX==datatype) { *indx = 19;
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
	    case '%':
		type = MOD;
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
		(void) sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("%s [%d]: syntax error in output",
			      __FILE__,__LINE__);
		break;
	    case POW:
	    case MOD:
	    case MULDIV:
	    case PLUSMIN:
		(void) sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("%s[%d]: syntax error in output",
			      __FILE__,__LINE__);
		(void) sf_pop(st1);
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
