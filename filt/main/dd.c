#include <stdio.h>

#include <rsf.h>

static size_t setfiledims (sf_file in, sf_file out);

static void ddbreak (sf_datatype itype, sf_datatype otype);

int main(int argc, char *argv[])
{
    size_t size, nin, nout, bufsiz=BUFSIZ;
    int line, ein, eout, n1, i, j, *ibuf;
    sf_file in, out;
    char *format, bufin[BUFSIZ], bufout[BUFSIZ];
    sf_datatype itype, otype;
    float *fbuf;
    float complex *cbuf;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output ("out");
    size = setfiledims(in,out);
    
    format = sf_getstring("data_format");
    if (NULL == format) sf_error("Specify data_format=");
	
    sf_setformat(out,format);

    if (!sf_histint(in,"esize",&ein)) ein=4;
    if (!sf_histint(out,"esize",&eout)) eout=4;

    itype = sf_gettype (in);
    otype = sf_gettype (out);
	
    if (SF_ASCII == sf_getform(out)) {
	if (!sf_getint("line",&line)) line=8;
	format = sf_getstring("format");
	sf_setaformat(format,line);
    }
    
    /* optimize buffer size */
    if (ein==0) {
	if (eout != 0) bufsiz /= eout;
    } else {
	bufsiz /= ein;
    }

    if (ein != 0 && eout != 0 && ein != eout && sf_histint(in,"n1",&n1)) {
	n1 = (n1*ein)/eout;
	sf_putint(out,"n1",n1);
    }

    while (size > 0) {
	nin = (bufsiz < size)? bufsiz:size;
	nout = (0==eout || 0==ein)? nin : nin*ein/eout;
	sf_read(bufin,ein,nin,in);
	if (itype != otype) {
	    switch (itype) {
		case SF_INT:
		    switch (otype) {
			case SF_FLOAT:
			    fbuf = (float*) bufout;
			    ibuf = (int*) bufin;
			    for (i=0, j=0; i < nin && j < nout; i++, j++) {
				*(fbuf++) = *(ibuf++); 
			    }
			    break;
			case SF_COMPLEX:
			    cbuf = (float complex*) bufout;
			    ibuf = (int*) bufin;
			    for (i=0, j=0; i < nin && j < nout; i++, j++) {
				*(cbuf++) = *(ibuf++); 
			    }
			    break;
			case SF_CHAR:
			case SF_INT:
			default:
			    ddbreak (itype,otype);
		    }
		    break;
		case SF_FLOAT:
		case SF_COMPLEX:
		case SF_CHAR:
		default:
		    ddbreak (itype,otype);
	    }
	    sf_write(bufout,eout,nout,out);
	} else {	    
	    sf_write(bufin,eout,nout,out);
	}
	size -= nin;
    }
    
    exit (0);
}

static void ddbreak (sf_datatype itype, sf_datatype otype)
{
    sf_error("Conversion from %d to %d"
	     "is unsupported",itype,otype);
}

static size_t setfiledims (sf_file in, sf_file out) 
{
    size_t size;
    int i, ni;
    char key[3];

    for (size=1, i=1; i <= SF_MAX_DIM; i++, size *= ni) {
	(void) snprintf(key,3,"n%d",i);
	if (!sf_getint(key,&ni)) {
	    if (!sf_histint(in,key,&ni)) break;
	} else {
	    sf_putint(out,key,ni);
	}
    }
    return size;
}
