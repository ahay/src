/* Create a header mask.

Mask is an integer data with ones and zeros. 
Ones correspond to header values between min and max.
The output can be used with sfheaderwindow.
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
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[]) {
    int *ibuf;
    size_t nsiz, nbuf, j;
    float *fbuf, min, max;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    sf_settype(out,SF_INT);

    nbuf = BUFSIZ/sizeof(float);
    fbuf = sf_floatalloc (nbuf);
    ibuf = sf_intalloc (nbuf);

    if (!sf_getfloat("min",&min)) min=-FLT_MAX;
    /* minimum header value */
    if (!sf_getfloat("max",&max)) max=+FLT_MAX;
    /* maximum header value */

    for (nsiz = sf_filesize (in); nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;

	sf_floatread(fbuf,nbuf,in);
	for (j=0; j < nbuf; j++) {
	    ibuf[j] = (fbuf[j] <= max && fbuf[j] >= min);
	}
	sf_intwrite(ibuf,nbuf,out);
    }

    exit(0);
}
	    
/* 	$Id: mask.c,v 1.7 2004/07/02 11:54:37 fomels Exp $	 */
