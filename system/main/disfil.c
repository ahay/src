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

static void head(int i, int cols, bool number);

int main (int argc, char* argv[])
{
    int cols, *ibuf, esize, i;
    short *sbuf;
    size_t bufsiz, nbuf, j;
    off_t size;
    const char* format;
    char *buf, *header, *trailer;
    float *fbuf;
    sf_complex *cbuf;
    unsigned char *ubuf;
    sf_file in;
    sf_datatype type;
    bool number;

    sf_init (argc,argv);
    
    in   = sf_input("in");
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
    header = sf_getstring("header");
    /* Optional header string to output before data */
    trailer = sf_getstring("trailer");
    /* Optional trailer string to output after data */

    esize = sf_esize(in);
    bufsiz = sf_bufsiz(in);
    if (0 != esize) bufsiz /= esize;
    size = sf_filesize(in);
    
    if (header != NULL) printf ("%s\n",header);   /* print header string */

    switch (type) {
	case SF_UCHAR:
	    if (0==cols) cols=10;
	    if (NULL==format) format = "%4d ";
	    ubuf = sf_ucharalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_ucharread (ubuf,nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    head(i,cols,number);
		    printf(format,ubuf[j]);
		}
	    }
	    printf("\n");
	    break;
	case SF_CHAR:
	    if (0==cols) cols=10;
	    if (NULL==format) format = "%4d ";
	    buf = sf_charalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_charread (buf,nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    head(i,cols,number);
		    printf(format,buf[j]);
		}
	    }
	    printf("\n");
	    break;
        case SF_SHORT:
	    if (0==cols) cols=10;
	    if (NULL==format) format = "%4d ";
	    sbuf = sf_shortalloc (bufsiz);
	    for (i=0; size > 0; size -= nbuf) {
		nbuf = (bufsiz < size)? bufsiz: size;
		sf_shortread (sbuf,nbuf,in);
		for (j=0; j < nbuf; j++, i++) {
		    head(i,cols,number);
		    printf(format,sbuf[j]);
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
		    head(i,cols,number);
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
		    head(i,cols,number);
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
		    head(i,cols,number);
		    printf(format,crealf(cbuf[j]),cimagf(cbuf[j]));
		}
	    }
	    printf("\n");
	    break;
	default:
	    sf_error("Unknown data type");
	    break;
    }

    if (trailer != NULL) printf ("%s\n",trailer);   /* print trailer string */


    exit (0);
}

static void head(int i, int cols, bool number) 
{
    if (0 == i%cols) {
	if (number) {
	    printf(i? "\n%4d: ":"%4d: ",i);
	} else if (i) {
	    printf("\n");
	}
    }
}

/* 	$Id: disfil.c 7791 2011-10-29 13:22:26Z sfomel $	 */
