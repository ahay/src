#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int resize, iesize, rsize, isize;
    size_t nleft, nbuf, i;
    sf_file real, imag, cmplx;
    char rbuf[BUFSIZ], ibuf[BUFSIZ], *cbuf, *rformat, *cformat;

    sf_init(argc,argv);
    if (3 != argc) 
	sf_error("wrong number of inputs: %d",argc-1);

    real = sf_input(argv[1]);
    imag = sf_input(argv[2]);
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
    sf_setformat(real,"raw");
    sf_setformat(imag,"raw");
    sf_setformat(cmplx,"raw");

    cbuf = sf_charalloc(2*BUFSIZ);

    for (nleft= (size_t) (rsize*resize); nleft > 0; nleft -= nbuf) {
	nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
	sf_read(rbuf,1,nbuf,real);
	sf_read(ibuf,1,nbuf,imag);
	for (i=0; i < nbuf; i += resize) {
	    memcpy(cbuf+2*i,       rbuf+i,(size_t) resize);
	    memcpy(cbuf+2*i+resize,ibuf+i,(size_t) resize);
	}
	sf_write(cbuf,1,2*nbuf,cmplx);
    }

    exit (0);
}
