#include <math.h>
#include <string.h>
#include <ctype.h>

#include <rsf.h>

static void check_compat (size_t nin, sf_file *in, int dim, const int *n);
static void parse (char* exp, int len, sf_file out,
		   sf_stack st1, sf_stack st2);
static void check (sf_stack st1, sf_stack st2);
static void evaluate (const sf_stack st, int nbuf, float** fbuf, float** fst);

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

enum {GRP, NUM, INDX, FUN, POW, MULDIV, PLUSMIN};

int main(int argc, char* argv[])
{
    int nin, i, n[SF_MAX_DIM], dim, nbuf, nsiz;
    size_t len;
    sf_file *in, out;
    sf_stack st1, st2;
    char *eq, *output, *key, *arg;
    float **fbuf, **fst;

    sf_init (argc,argv);
    in = (sf_file*) sf_alloc ((size_t) argc-1,sizeof(sf_file));
    out = sf_output ("out");

    nin=0;
    for (i=1; i< argc; i++) { /* collect inputs */
	arg = argv[i];
	eq =  strchr(arg,'=');
	if (NULL == eq) continue; /* not a parameter */
	if (0 == strncmp(arg,"output",6)) continue; /* not a file */

	len = (size_t) (eq-arg);
	key = sf_charalloc(len+1);
	strncpy(key,arg,len);
	key[len]='\0';
	 
	in[nin] = sf_input(key);
	sf_putint(out,key,nin);
	nin++;
	free(key);
    }

    dim = sf_filedims(in[0],n);
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }
    
    check_compat(nin,in,dim,n);

    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");

    len = strlen(output);
    st1 = sf_stack_init (len);
    st2 = sf_stack_init (len);

    /* algorithm is in st2 */
    parse (output,len,out,st1,st2);
    len = sf_stack_get(st2);

    check (st1,st2);
    sf_stack_close (st1);

    nbuf = BUFSIZ/sizeof(float);

    fbuf = sf_floatalloc2(nbuf,nin);
    fst  = sf_floatalloc2(nbuf,len+2);
    
    sf_fileflush(out,in[0]);

    for (; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;
	for (i=0; i < nin; i++) {
	    sf_read(fbuf[i],sizeof(float),nbuf,in[i]);
	}

	sf_stack_set (st2,len);
	evaluate (st2, nbuf, fbuf, fst);

	sf_write(fst[1],sizeof(float),nbuf,out);
    }
    
    exit(0);
}

static void evaluate (const sf_stack st, int nbuf, float** fbuf, float** fst)
{
    char *op;
    int *indx, i;
    float *num, f, *farr;
    func fun;

    while (sf_full (st)) {
	switch (sf_top (st)) {
	    case NUM: 
		farr = *(++fst);
		num = (float*) sf_pop(st);
		f = *num;
		for (i=0; i < nbuf; i++) { farr[i] = f; }
		break;
	    case INDX:
		farr = *(++fst);
		indx = (int*) sf_pop(st);
		/* convert index to number */
		num = *(fbuf + (*indx));
		for (i=0; i < nbuf; i++) { farr[i] = num[i]; }
		break;
	    case FUN:
		indx = (int*) sf_pop(st);
		fun = functable[*indx];
		farr = *fst;
		for (i=0; i < nbuf; i++) { farr[i] = fun(farr[i]); }
		break;
	    case POW:
		sf_pop(st);
		num = *fst;
		farr = *(--fst);
		for (i=0; i < nbuf; i++) { farr[i] = powf(farr[i],num[i]); }
		break;
	    case MULDIV:
		op = (char*) sf_pop(st);
		num = *fst;
		farr = *(--fst);
		if ('*' == *op) {
		    for (i=0; i < nbuf; i++) { farr[i] *= num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { farr[i] /= num[i]; }
		}
		break;
	    case PLUSMIN:
		op = (char*) sf_pop(st);
		num = *fst;
		farr = *(--fst);
		if ('+' == *op) {
		    for (i=0; i < nbuf; i++) { farr[i] += num[i]; }
		} else {
		    for (i=0; i < nbuf; i++) { farr[i] -= num[i]; }
		}
		break;
	    default:
		sf_error ("syntax error in output");
		break;
	}
    }
}

static void parse (char* output, int len, sf_file out, 
		   sf_stack st1, sf_stack st2)
{
    int i, j, keylen, *indx, type=-1, top;
    char *key, c, c2;
    float *num;

    for (i=0; i < len; i++) {
	c = output[i];
	
	if (isspace(c)) continue; /* skip white space */ 

	/* handle parentheses */

	if ('(' == c) {
	    sf_push(st2,&c,GRP);
	    continue;
	}

	if (')' == c) {
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
		sf_error("unbalanced ')', position %d in output",i);
	    continue;
	}

	if ('.' == c || isdigit(c)) { /* number */
	    for (j=i+1; j < len; j++) {
		c2 = output[j];
		if ('.' != c2 || !isdigit(c2)) break;
	    }
	    keylen = j-i;
	    key = sf_charalloc(keylen+1);
	    strncpy(key,output+i,keylen);
	    key[keylen]='\0';
 
	    num = sf_floatalloc(1);
	    if (0 == sscanf(key,"%f",num)) 
		sf_error("ill-formated number: %s, position %d in output",key,i);
	    free(key);

	    sf_push(st1,num,NUM); /* st1 is out */

	    i=j-1;
	    continue;
	}

	if (isalpha (c)) { /* identifier */
	    for (j=i+1; j < len; j++) {
		if (!isalnum(output[j])) break;
	    }

	    keylen = j-i;
	    key = sf_charalloc(keylen+1);
	    strncpy(key,output+i,keylen);
	    key[keylen]='\0';

	    indx = sf_intalloc(1);

	    if (sf_histint(out,key,indx)) {
		sf_push(st1,indx,INDX); 
	    } else {
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
		    sf_error("unrecognized identifier: "
			     "%s, position %d in output",key,i);
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
		type = PLUSMIN;
		break;
	    case '*':
	    case '/':
		type = MULDIV;
		break;
	    case '^':
		type = POW;
		break;
	    default:
		sf_error("wrong character: %c, position %d in output",c);
		break;
	}

	while (sf_full(st2)) {
	    top = sf_top (st2);
	    if (GRP==top || top > type) break; /* compare precedence */
	    sf_push(st1,sf_pop(st2),top);
	}
	
	sf_push(st2,output+i,type);
    }
	
    /* push operators into output (reverse polish=) */
    while (sf_full(st2)) {
	top = sf_top (st2);
	if (GRP == top)
	    sf_error("unbalanced '(' in output");
	sf_push(st1,sf_pop(st2),top);
    }

    /* arrange in reverse order for computation */
    while (sf_full(st1)) {
	top = sf_top(st1);
	sf_push(st2,sf_pop(st1),top);
    }
}

static void check (sf_stack st1, sf_stack st2)
{
    int top;

    while (sf_full (st2)) {
	switch (sf_top (st2)) {
	    case NUM: 
	    case INDX:
		sf_push(st1,sf_pop(st2),NUM);
		break;
	    case FUN:
		sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("syntax error in output");
		break;
	    case POW:
	    case MULDIV:
	    case PLUSMIN:
		sf_pop(st2);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("syntax error in output");
		sf_pop(st1);
		top = sf_top (st1);
		if (NUM != top && INDX != top) 
		    sf_error ("syntax error in output");
		break;
	    default:
		sf_error ("syntax error in output");
		break;
	}
    }
    top = sf_top (st1);
    if (0 != sf_stack_get(st1) || (NUM != top && INDX != top)) 
	sf_error ("syntax error in output");
}

static void check_compat (size_t nin, sf_file *in, int dim, const int *n) 
{
    int ni, id;
    size_t i;
    float f, fi;
    char key[3];
    const float tol=1.e-5;
    
    if (SF_FLOAT != sf_gettype(in[0])) 
	sf_error("Need float input");
    for (i=1; i < nin; i++) {
	if (SF_FLOAT != sf_gettype(in[i])) 
	    sf_error("Need float input");
	for (id=1; id <= dim; id++) {
	    (void) snprintf(key,3,"n%d",id);
	    if (!sf_histint(in[i],key,&ni) || ni != n[id-1])
		sf_error("%s mismatch: need %d",key,n[id-1]);
	    (void) snprintf(key,3,"d%d",id);
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || (fabs(fi-f) > tol*fabs(f))))
		sf_warning("%s mismatch: need %g",key,f);
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || (fabs(fi-f) > tol*fabs(f))))
		sf_warning("%s mismatch: need %g",key,f);
	}
    }
}
