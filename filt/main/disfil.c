#include <rsf.h>

int main (int argc, char* argv[])
{
    int cols, *ibuf, esize, i;
    size_t bufsiz = BUFSIZ, nbuf, size, j;
    char* format, *buf;
    float *fbuf;
    float complex *cbuf;
    sf_file in;
    sf_datatype type;
    bool number;

    sf_init (argc,argv);
    
    in = sf_input ("in");
    type = sf_gettype(in);

    if (!sf_getbool("number",&number)) number=true;
    if (!sf_getint("col",&cols)) cols=0; 
    format = sf_getstring("format");

    if (!sf_histint(in,"esize",&esize)) esize=4;
    if (0 != esize) bufsiz /= esize;
    size = (size_t) sf_filesize(in);

    switch (type) {
	case SF_CHAR:
	    if (0==cols) cols=10;
	    if (NULL==format) format = "%4d ";
	    buf = sf_charalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_read (buf,sizeof(char),nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    if (number && 0 == i%cols) printf(i? "\n%4d: ":"%4d: ",i);
		    printf(format,buf[j]);
		}
	    }
	    printf("\n");
	    break;
	case SF_INT:
	    if (0==cols) cols=10;
	    if (NULL==format) format = "%4d ";
	    ibuf = sf_intalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_read (ibuf,sizeof(int),nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    if (number && 0 == i%cols) printf(i? "\n%4d: ":"%4d: ",i);
		    printf(format,ibuf[j]);
		}
	    }
	    printf("\n");
	    break;
	case SF_FLOAT:
	    if (0==cols) cols=5;
	    if (NULL==format) format = "%13.4g";
	    fbuf = sf_floatalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_read (fbuf,sizeof(float),nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    if (0 == i%cols) printf(i? "\n%4d: ":"%4d: ",i);
		    printf(format,fbuf[j]);
		}
	    }
	    printf("\n");
	    break;
	case SF_COMPLEX:
	    if (0==cols) cols=3;
	    if (NULL==format) format = "%10.4g,%10.4gi";
	    cbuf = sf_complexalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_read (cbuf,sizeof(float complex),nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    if (0 == i%cols) printf(i? "\n%4d: ":"%4d: ",i);
		    printf(format,crealf(cbuf[j]),cimagf(cbuf[j]));
		}
	    }
	    printf("\n");
	    break;
	default:
	    sf_error("Unknown data type");
	    break;
    }

    exit (0);
}
