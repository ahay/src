/* Print out data values.

Alternatively, use sfdd and convert to ASCII form.
*/
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
    /* If number the elements */
    if (!sf_getint("col",&cols)) cols=0;
    /* Number of columns.
       The default depends on the data type:
       10 for int and char,
       5 for float,
       3 for complex */
    format = sf_getstring("format");
    /* Format for numbers (printf-style).
       The default depends on the data type:
       "%4d " for int and char,
       "%13.4g" for float,
       "%10.4g,%10.4gi" for complex */

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
		sf_charread (buf,nbuf,in);
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
		sf_intread (ibuf,nbuf,in);
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
		sf_floatread (fbuf,nbuf,in);
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
		sf_complexread (cbuf,nbuf,in);
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

/* 	$Id: disfil.c,v 1.6 2004/07/02 11:54:37 fomels Exp $	 */
