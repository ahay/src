/* Convert real data to complex (by adding zero imaginary part).

Takes: < real.rsf > cmplx.rsf

See also: sfcmplx
*/

#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int esize, size;
    size_t i, nleft, nbuf, e_size=4;
    sf_file real, cmplx;
    char rbuf[BUFSIZ], *cbuf, *rformat, *cformat;

    sf_init(argc,argv);
    real = sf_input("in");
    cmplx = sf_output ("out");

    if (SF_FLOAT != sf_gettype(real)) sf_error("wrong input type");
    
    rformat = sf_histstring(real,"data_format");
    cformat = sf_charalloc(strlen(rformat)+1
			   -strlen("float")
			   +strlen("complex"));
    strncpy(cformat,rformat,strlen(rformat));
    strcpy(strstr(cformat,"float"),"complex");
    sf_setformat(cmplx,cformat);

    if (sf_histint(real,"esize",&esize)) {
	if (esize <= 0) {
	    sf_error("wrong esize=%d",esize);
	} else {
	    e_size=(size_t) esize;
	}
    }

    size = sf_filesize (real);

    sf_fileflush(cmplx,real);
    sf_setformat(real,"raw");
    sf_setformat(cmplx,"raw");

    cbuf = sf_charalloc(2*BUFSIZ);
    for (i=0; i < BUFSIZ; i += e_size) {
	memset(cbuf+2*i+e_size,0,e_size);
    }

    for (nleft=size*e_size; nleft > 0; nleft -= nbuf) {
	nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	sf_read(rbuf,1,nbuf,real);
	for (i=0; i < nbuf; i += e_size) {
	    memcpy(cbuf+2*i,rbuf+i,e_size);
	}
	sf_write(cbuf,1,2*nbuf,cmplx);
    }

    exit (0);
}

/* 	$Id: rtoc.c,v 1.2 2003/09/29 14:34:56 fomels Exp $	 */
