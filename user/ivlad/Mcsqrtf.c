/* Complex square root. Good example of I/O loop for applying a function.
Realized after I wrote this program that the sqrt function in sfmath does the
same job, but keeping it around as a simple example of buffer I/O. */
/*
  Copyright (C) 2012 Ioan Vlad

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

int main (int argc, char *argv[]) {

    sf_file in=NULL, out=NULL;
    sf_complex* buf=NULL;
    off_t nsiz;
    size_t bufsiz, ibuf, nbuf, nleft;    

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in))
        sf_error( "Need complex input" );

    nsiz = sf_leftsize(in, 0);
    bufsiz = sf_bufsiz(in) / sf_esize(in);
    buf = sf_complexalloc(bufsiz);

    for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
	    nbuf = (bufsiz < nleft)? bufsiz: nleft;
        sf_complexread(buf, nbuf, in);
        for (ibuf=0; ibuf < nbuf; ibuf++) {
            buf[ibuf] = csqrtf(buf[ibuf]);
        }
        sf_complexwrite(buf, nbuf, out);
    }

    exit(0);
}
