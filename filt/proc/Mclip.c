/* Clip the data.
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

int main(int argc, char* argv[])
{
    int i, n, nbuf;
    float clip, *trace;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    n = sf_filesize(in);

    /* parameter from the command line (i.e. clip=1.5 ) */
    if (!sf_getfloat("clip",&clip)) sf_error("Need clip=");

    /* allocate floating point array */
    nbuf = BUFSIZ/sizeof(float);
    trace = sf_floatalloc (nbuf);

    /* loop over traces */
    for (n = sf_filesize(in); n > 0; n -= nbuf) {
	if (nbuf > n) nbuf=n;

	/*read a trace */
	sf_floatread(trace,nbuf,in);

	/* loop over samples */
	for (i=0; i < n; i++) {
	    if      (trace[i] >  clip) trace[i]= clip;
	    else if (trace[i] < -clip) trace[i]=-clip;
	}
    
	/* write a trace */
	sf_floatwrite(trace,nbuf,out);
    }
    
    sf_close();
    exit(0);
}

/* 	$Id: Mclip.c,v 1.2 2004/06/25 08:41:19 fomels Exp $	 */

#ifndef lint
static char vcid[] = "$Id: Mclip.c,v 1.2 2004/06/25 08:41:19 fomels Exp $";
#endif /* lint */
