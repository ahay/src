/* Create a complex dataset from its real and imaginary parts.
   
Takes: real.rsf imag.rsf

There has to be only two input files specified and no additional parameters.
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

#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int resize, iesize;
    off_t nleft, nbuf, i, rsize, isize;
    sf_file real=NULL;   /* input */
    sf_file imag=NULL;   /* input */
    sf_file cmplx=NULL; /* output */
    char rbuf[BUFSIZ], ibuf[BUFSIZ], *cbuf;

    sf_init(argc,argv);

    /* the first two non-parameters are real and imag files */
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
    if (NULL == imag) {
	if (NULL == real) sf_error ("not enough input");
	/* if only one input, real is in stdin */
	imag = real;
	real = sf_input("in");
    }
    cmplx = sf_output ("out");

    if (SF_FLOAT != sf_gettype(real) ||
	SF_FLOAT != sf_gettype(imag))
	sf_error("wrong input type");

    sf_settype(cmplx,SF_COMPLEX);

    resize = sf_esize(real);
    iesize = sf_esize(imag);
    if (resize != iesize) sf_error("esize mismatch: %d != %d",resize,iesize);

    rsize = sf_filesize (real);
    isize = sf_filesize (imag);
    if (rsize != isize) sf_error("size mismatch: %d != %d",rsize,isize);

    sf_fileflush(cmplx,real);
    sf_setform( real,SF_NATIVE);
    sf_setform( imag,SF_NATIVE);
    sf_setform(cmplx,SF_NATIVE);

    cbuf = sf_charalloc(2*BUFSIZ);

    for (nleft= (size_t) (rsize*resize); 
	 nleft > 0; nleft -= nbuf) {
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
