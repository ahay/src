/* Convert between different formats. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <stdio.h>

#include <rsf.h>

static size_t setfiledims (sf_file in, sf_file out);

static void ddbreak (sf_datatype itype, sf_datatype otype);

int main(int argc, char *argv[])
{
    size_t size, nin, nout, bufsiz=BUFSIZ;
    int line, ein, eout, n1, i, j, *ibuf;
    sf_file in, out;
    char *form, *type, *format, bufin[BUFSIZ], bufout[BUFSIZ];
    sf_datatype itype, otype;
    float *fbuf;
    float complex *cbuf;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output ("out");
    size = setfiledims(in,out);
    
    form = sf_getstring("form");
    /* ascii, native, xdr */
    type = sf_getstring("type");
    /* int, float, complex */
    if (NULL == form && NULL == type) sf_error("Specify form= or type=");
	
    if (NULL != form) {
	switch (form[0]) {
	    case 'a':
		sf_setform(out,SF_ASCII);
		if (!sf_getint("line",&line)) line=8;
		/* Number of numbers per line (for conversion to ASCII) */
		format = sf_getstring("format");
		/* Element format (for conversion to ASCII) */ 
		sf_setaformat(format,line);
		break;
	    case 'n':
		sf_setform(out,SF_NATIVE);
		break;
	    case 'x':
		sf_setform(out,SF_XDR);
		break;
	    default:
		sf_error("Unsupported form=\"%s\"",form);
		break;
	}
    }

    itype = sf_gettype (in);
    otype = itype;
    if (NULL != type) {
	switch (type[0]) {
	    case 'i':
		otype = SF_INT;
		break;
	    case 'f':
		otype = SF_FLOAT;
		break;
	    case 'c':
		otype = SF_COMPLEX;
		break;
	    default:
		sf_error("Unsupported type=\"%s\"",type);
		break;
	}
	sf_settype(out,otype);
    } 
	
    if (!sf_histint(in,"esize",&ein)) ein=4;
    if (!sf_histint(out,"esize",&eout)) eout=4;
    
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
	switch (itype) {
	    case SF_INT:
		ibuf = (int*) bufin;
		sf_intread(ibuf,nin,in);
		switch (otype) {
		    case SF_INT:
			sf_intwrite(ibuf,nout,out);
			break;
		    case SF_FLOAT:
			fbuf = (float*) bufout;
			for (i=j=0; i < nin && j < nout; i++, j++) {
			    fbuf[j] = ibuf[i]; 
			}
			sf_floatwrite(fbuf,nout,out);
			break;
		    case SF_COMPLEX:
			cbuf = (float complex*) bufout;
			for (i=j=0; i < nin && j < nout; i+=2, j++) {
			    cbuf[j] = ibuf[i]+I*ibuf[i+1]; 
			}
			sf_complexwrite(cbuf,nout,out);
			break;
		    case SF_CHAR:
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_FLOAT:
		fbuf = (float*) bufin;
		sf_floatread(fbuf,nin,in);
		switch (otype) {
		    case SF_FLOAT:
			sf_floatwrite(fbuf,nout,out);
			break;
		    case SF_INT:
			ibuf = (int*) bufout;
			for (i=j=0; i < nin && j < nout; i++, j++) {
			    ibuf[i] = fbuf[j]; 
			}
			sf_intwrite(ibuf,nout,out);
			break;
		    case SF_COMPLEX:
			cbuf = (float complex*) bufout;
			for (i=j=0; i < nin && j < nout; i+=2, j++) {
			    cbuf[j] = fbuf[i] + I*fbuf[i+1];
			}
			sf_complexwrite(cbuf,nout,out);
			break;			    
		    case SF_CHAR:
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_COMPLEX:
	    case SF_CHAR:
	    default:
		ddbreak (itype,otype);
		break;
	}
	size -= nin;
    }
    
    exit (0);
}

static void ddbreak (sf_datatype itype, sf_datatype otype)
{
    const char* types[]={"char","int","float","complex"};

    sf_error("Conversion from %s to %s"
	     " is unsupported",types[itype],types[otype]);
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

/* 	$Id$	 */
