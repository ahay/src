/* Combine several datasets by interleaving.

Takes: [< file0.rsf] file1.rsf file2.rsf ...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat (size_t esize, 
			  int nin, sf_file* in, int dim, const off_t *n);

int main (int argc, char* argv[])
{
    int i, nin, axis, dim;
    off_t n[SF_MAX_DIM], n1, i2, n2, nleft;
    size_t nbuf, esize;
    sf_file *in=NULL, out=NULL;
    char key[3], buf[BUFSIZ];

    sf_init(argc,argv);
    in = (sf_file*) sf_alloc ((size_t) argc,sizeof(sf_file));

    if (!sf_stdin()) { /* no input file in stdin */
	nin=0;
    } else {
	in[0] = sf_input("in");
	nin=1;
    }

    for (i=1; i< argc; i++) { /* find other */
	if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	in[nin] = sf_input(argv[i]);
	nin++;
    }
    if (0==nin) sf_error ("no input");

    dim = sf_largefiledims(in[0],n);
    if (!sf_getint("axis",&axis)) axis=3;
    /* Axis for interleaving */
    if (axis > dim && axis <=0) 
	sf_error("axis=%d is not in the range [1,%d]",axis,dim);

    esize=sf_esize(in[0]);
    check_compat(esize,nin,in,dim,n);

    n1 = esize;
    for (i=1; i < axis; i++) {
	n1 *= n[i-1];
    }
    for (n2=1; i <= dim; i++) {
	n2 *= n[i-1];
    }

    out = sf_output ("out");
    sprintf(key,"n%d",axis);
    sf_putint(out,key,n[axis-1]*nin);
    sf_setformat(out,sf_histstring(in[0],"data_format"));

    sf_fileflush(out,in[0]);
    for (i=0; i < nin; i++) {
	sf_setform(in[i],SF_NATIVE);
    }
    sf_setform(out,SF_NATIVE);

    for (i2=0; i2 < n2; i2++) {
	for (i=0; i < nin; i++) {
	    for (nleft=n1; nleft > 0; nleft -= nbuf) {
		nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
		sf_charread  (buf,nbuf,in[i]);
		sf_charwrite (buf,nbuf,out);
	    }
	}
    }


    exit(0);
}

static void check_compat (size_t esize, 
			  int nin, sf_file *in, int dim, const off_t *n)
{
    int ni, id, i;
    char key[3];

    for (i=1; i < nin; i++) {
	if (esize != sf_esize(in[i]))
	    sf_error ("esize mismatch: need %d",esize);
	for (id=0; id < dim; id++) {
	    snprintf(key,3,"n%d",id+1);
	    if (!sf_histint(in[i],key,&ni) || ni != n[id])
		sf_error("%s mismatch: need %d",key,n[id]);
	}
    }
}
