/* Combine two datasets by interleaving.

Alternatively, specify: in.rsf other.rsf > out.rsf
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat (sf_file other, int esize, int dim, const int *n);

int main (int argc, char* argv[])
{
    int i, axis, n1, n2, i2, esize, dim, n[SF_MAX_DIM];
    size_t nbuf, nleft;
    sf_file in, other, out;
    char key[3], buf[BUFSIZ];
    
    sf_init(argc,argv);

    if (0 != isatty(fileno(stdin))) { /* no input file in stdin */
	in = NULL;
	for (i=1; i< argc; i++) { /* find in */
	    if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	    in = sf_input(argv[i]);
	    break;
	}
	if (NULL == in) sf_error ("no input");
    } else {
	in = sf_input("in");
	i=0;
    }

    if (NULL != sf_getstring("other")) {
	/* Other dataset */
	other = sf_input("other");
    } else {
	other = NULL;
	for (i++; i< argc; i++) { /* find other */
	    if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	    other = sf_input(argv[i]);
	    break;
	}
	if (NULL == other) sf_error ("not enough input");
    }
   
    dim = sf_filedims(in,n);
    if (!sf_getint("axis",&axis)) axis=3;
    /* Axis for interleaving */
    if (axis > dim && axis <=0) 
	sf_error("axis=%d is not in the range [1,%d]",axis,dim);

    if (!sf_histint(in,"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("wrong esize=%d",esize);
    }
    check_compat(other,esize,dim,n);
    
    n1 = esize;
    for (i=1; i < axis; i++) {
	n1 *= n[i-1];
    }
    for (n2=1; i <= dim; i++) {
	n2 *= n[i-1];
    }
    
    sprintf(key,"n%d",axis);

    out = sf_output ("out");
    sf_putint(out,key,n[axis-1]*2);
    sf_setformat(out,sf_histstring(in,"data_format"));

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(other,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    for (i2=0; i2 < n2; i2++) {
	for (nleft=(size_t) n1; nleft > 0; nleft -= nbuf) {
	    nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	    sf_charread  (buf,nbuf,in);
	    sf_charwrite (buf,nbuf,out);
	}
	for (nleft=(size_t) n1; nleft > 0; nleft -= nbuf) {
	    nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	    sf_charread  (buf,nbuf,other);
	    sf_charwrite (buf,nbuf,out);
	}
    }
    
    exit(0);
}

static void check_compat (sf_file other, int esize, int dim, const int *n)
{
    int ni, id;
    char key[3];
    
    if (!sf_histint(other,"esize",&ni) || ni != esize)
	sf_error ("esize mismatch: need %d",esize);
    for (id=0; id < dim; id++) {
	snprintf(key,3,"n%d",id+1);
	if (!sf_histint(other,key,&ni) || ni != n[id])
	    sf_error("%s mismatch: need %d",key,n[id]);
    }
}

/* 	$Id: interleave.c,v 1.8 2004/07/02 11:54:37 fomels Exp $	 */
