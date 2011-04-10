/* Extract real (sfreal) or imaginary (sfimag) part of a complex dataset. */
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
    int esize, shift;
    off_t size;
    size_t i, nleft, nbuf, e_size, len, bufsiz;
    sf_file real, cmplx;
    char *rbuf, *cbuf, *rformat, *cformat, *prog;

    sf_init(argc,argv);
    cmplx = sf_input ( "in");
    real  = sf_output("out");

    if (SF_COMPLEX != sf_gettype(cmplx))
      {
	
	sf_error("wrong input type");
      }
    
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
    sf_setform(real ,SF_NATIVE);
    sf_setform(cmplx,SF_NATIVE);
    
    prog = sf_getprog();
    if (       NULL != strstr(prog,"real")) {
	shift=0;
    } else if (NULL != strstr(prog,"imag")) {
	shift=esize;
    } else {
	sf_warning("neither real nor imag, assume real");
	shift=0;
    }

    bufsiz = sf_bufsiz(cmplx);
    rbuf = sf_charalloc(bufsiz);
    cbuf = sf_charalloc(2*bufsiz);
    
    for (nleft=size*e_size; nleft > 0; nleft -= nbuf) {
	nbuf = (bufsiz < nleft)? bufsiz: nleft;
	sf_charread(cbuf,2*nbuf,cmplx);
	for (i=0; i < nbuf; i += e_size) {
	    memcpy(rbuf+i,cbuf+2*i+shift,e_size);
	}
	sf_charwrite(rbuf,nbuf,real);
    }
    

    exit (0);
}
