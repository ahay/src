/* Convert real data to complex (by adding zero imaginary part).

See also: sfcmplx
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
    int size;
    off_t i, nleft, nbuf, e_size;
    sf_file real, cmplx;
    bool pair = false;
    char rbuf[BUFSIZ], *cbuf, *rformat, *cformat;

    sf_init(argc,argv);
    real = sf_input("in");
    cmplx = sf_output ("out");

    if (SF_FLOAT != sf_gettype(real)) sf_error("wrong input type");

    if (!sf_getbool("pair",&pair)) pair = false;
    /* y - use odd elements for real part and even ones for imaginary part */
    if (pair) {
        if (!sf_histint(real,"n1",&size)) sf_error("No n1= in input");
        if ((size % 2)) sf_error("n1= in input should be even");
        size /= 2;
        if (size != 1)
            sf_putint(cmplx,"n1",size);
        else
            sf_unshiftdim(real,cmplx,1);
    }

    rformat = sf_histstring(real,"data_format");
    cformat = sf_charalloc(strlen(rformat)+1
                           -strlen("float")
                           +strlen("complex"));
    strncpy(cformat,rformat,strlen(rformat));
    strcpy(strstr(cformat,"float"),"complex");
    sf_setformat(cmplx,cformat);

    e_size= sf_esize(real);
    size = sf_filesize (real);

    sf_fileflush(cmplx,real);
    sf_setform(real,SF_NATIVE);
    sf_setform(cmplx,SF_NATIVE);

    cbuf = sf_charalloc(2*BUFSIZ);
    for (i=0; i < BUFSIZ; i += e_size) {
        memset(cbuf+2*i+e_size,0,e_size);
    }

    for (nleft=size*e_size; nleft > 0; nleft -= nbuf) {
        nbuf = (BUFSIZ < nleft)? BUFSIZ: nleft;
        sf_charread(rbuf,nbuf,real);
        if (pair) {
            sf_charwrite(rbuf,nbuf,cmplx);
        } else {
            for (i=0; i < nbuf; i += e_size) {
                memcpy(cbuf+2*i,rbuf+i,e_size);
            }
            sf_charwrite(cbuf,2*nbuf,cmplx);
        }
    }

    exit (0);
}

/*         $Id$         */
