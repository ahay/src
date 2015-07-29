/* Computes number of dimensions and their values
Wrapper for sf_filedims. */
/*
  Copyright (C) 2007 Ioan Vlad

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

int
main (int argc, char *argv[])
{
    int i;
    int ndims;
    off_t nlarge[SF_MAX_DIM];
    int n[SF_MAX_DIM];
    bool parform, large;
    sf_file in=NULL;

    sf_init(argc,argv);

    in = sf_input("in");

    if (!sf_getbool("large", &large))
        large = false;
    /* if y, file with large dimensions. */

    if (large) {
	ndims = sf_largefiledims(in, nlarge);
    } else {
	ndims = sf_filedims(in, n);
    }

    if (!sf_getbool("parform", &parform))
        parform = true;
    /* If y, print out parameter=value. If n, print out value. */

    if (parform) printf( "ndims=" );

    printf( "%d", ndims );

    if (parform)
        printf( "\nn=" );
    else
        printf( ":" );

    for (i=0; i<ndims; i++) {
	if (large) {
#if defined(__cplusplus) || defined(c_plusplus)
	    printf( "%ld", (long) n[i] );
#else
	    printf( "%lld", (long long) n[i] );
#endif
	} else {
	    printf( "%d", n[i] );
	}
        if (i<ndims-1) printf(",");
    }

    printf( "\n" );


    exit(0);
}
