/* Extract real (sfreal) or imaginary (sfimag) part of a complex dataset.

Takes: < cmplx.rsf > real.rsf
*/

#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int esize, size, shift;
    size_t i, nleft, nbuf, e_size, len;
    sf_file real, cmplx;
    char rbuf[BUFSIZ], *cbuf, *rformat, *cformat, *prog;

    sf_init(argc,argv);
    cmplx = sf_input("in");
    real = sf_output ("out");

    if (SF_COMPLEX != sf_gettype(cmplx)) sf_error("wrong input type");
    
    cformat = sf_histstring(cmplx,"data_format");
    len = strlen(cformat)+1;
    rformat = sf_charalloc(len);
    memcpy(rformat,cformat,len);
    strcpy(strstr(rformat,"complex"),"float");
    sf_setformat(real,rformat);

    if (!sf_histint(real,"esize",&esize)) esize=4;
    if (esize <= 0) sf_error("wrong esize=%d",esize);
    e_size = (size_t) esize;

    size = sf_filesize (cmplx);
    
    sf_fileflush(real,cmplx);
    sf_setformat(real,"raw");
    sf_setformat(cmplx,"raw");

    prog = sf_getprog();
    if (NULL != strstr(prog,"real")) {
	shift=0;
    } else if (NULL != strstr(prog,"imag")) {
	shift=esize;
    } else {
	sf_warning("neither real nor imag, assume real");
	shift=0;
    }

    cbuf = sf_charalloc(2*BUFSIZ);
    
    for (nleft=size*e_size; nleft > 0; nleft -= nbuf) {
	nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	sf_read(cbuf,1,2*nbuf,cmplx);
	for (i=0; i < nbuf; i += e_size) {
	    memcpy(rbuf+i,cbuf+2*i+shift,e_size);
	}
	sf_write(rbuf,1,nbuf,real);
    }
    sf_fileclose(real);

    exit (0);
}

/* 	$Id: real.c,v 1.2 2003/09/29 14:34:56 fomels Exp $	 */
