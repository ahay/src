/* Create a complex dataset from its real and imaginary parts.

Takes: real.rsf imag.rsf > cmplx.rsf

There has to be only two input files specified and no additional parameters.
*/

#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int resize, iesize, rsize, isize;
    size_t nleft, nbuf, i;
    sf_file real=NULL, imag=NULL, cmplx;
    char rbuf[BUFSIZ], ibuf[BUFSIZ], *cbuf, *rformat, *cformat;

    sf_init(argc,argv);

    /* the first two non-parameters are real and imaginary files */
    for (i=1; i< argc; i++) { 
	if (NULL == strchr(argv[i],'=')) {
	    if (NULL == real) {
		real = sf_input (argv[i]);
	    } else {
		imag = sf_input (argv[i]);
		break;
	    }
	}
    }
    if (NULL == real || NULL == imag)
	sf_error ("not enough input");
    cmplx = sf_output ("out");

    if (SF_FLOAT != sf_gettype(real) ||
	SF_FLOAT != sf_gettype(imag))
	sf_error("wrong input type");
    
    rformat = sf_histstring(real,"data_format");
    cformat = sf_charalloc(strlen(rformat)+1
			   -strlen("float")
			   +strlen("complex"));
    memcpy(cformat,rformat,strlen(rformat));
    strcpy(strstr(cformat,"float"),"complex");
    sf_setformat(cmplx,cformat);

    if (!sf_histint(real,"esize",&resize)) resize=4;
    if (!sf_histint(imag,"esize",&iesize)) iesize=4;
    if (resize != iesize) sf_error("esize mismatch: %d != %d",resize,iesize);
    if (resize <= 0) sf_error("wrong esize=%d",resize);
    
    rsize = sf_filesize (real);
    isize = sf_filesize (imag);
    if (rsize != isize) sf_error("size mismatch: %d != %d",rsize,isize);

    sf_fileflush(cmplx,real);
    sf_setform(real,SF_NATIVE);
    sf_setform(imag,SF_NATIVE);
    sf_setform(cmplx,SF_NATIVE);

    cbuf = sf_charalloc(2*BUFSIZ);

    for (nleft= (size_t) (rsize*resize); nleft > 0; nleft -= nbuf) {
	nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	sf_charread(rbuf,nbuf,real);
	sf_charread(ibuf,nbuf,imag);
	for (i=0; i < nbuf; i += resize) {
	    memcpy(cbuf+2*i,       rbuf+i,(size_t) resize);
	    memcpy(cbuf+2*i+resize,ibuf+i,(size_t) resize);
	}
	sf_charwrite(cbuf,2*nbuf,cmplx);
    }

    exit (0);
}

/* 	$Id: cmplx.c,v 1.4 2004/04/19 21:51:35 fomels Exp $	 */


