/* Combine two datasets by interleaving.

Takes: [< file.rsf | file.rsf] [other=other.rsf | other.rsf] > out.rsf
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
    sf_setformat(in,"raw");
    sf_setformat(other,"raw");
    sf_setformat(out,"raw");

    for (i2=0; i2 < n2; i2++) {
	for (nleft=(size_t) n1; nleft > 0; nleft -= nbuf) {
	    nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	    sf_read  (buf,1,nbuf,in);
	    sf_write (buf,1,nbuf,out);
	}
	for (nleft=(size_t) n1; nleft > 0; nleft -= nbuf) {
	    nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	    sf_read  (buf,1,nbuf,other);
	    sf_write (buf,1,nbuf,out);
	}
    }
    
    sf_close();
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

/* 	$Id: interleave.c,v 1.5 2004/03/22 05:43:24 fomels Exp $	 */
