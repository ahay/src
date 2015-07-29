/* dumps header information to the standard output.
   Extended sffiledims.*/
/*
  Copyright (C) 2013 Henryk Modzelewski

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
    bool parform, large, all;
    char key[16];
    float fdummy;
    sf_file in=NULL;

    sf_init(argc,argv);

    in = sf_input("in");

    if (!sf_getbool("large", &large))
        large = false;
    /* if y, file with large dimensions. */

    if (!sf_getbool("parform", &parform))
        parform = true;
    /* If y, print out parameter=value. If n, print out value. */

    if (!sf_getbool("all", &all))
        all = false;
    /* If y, print all values, icluding singleton dimensions.
       If n, drop trailing singleteon dimensions.*/

    if (large) {
	ndims = sf_largefiledims(in, nlarge);
    } else {
	ndims = sf_filedims(in, n);
    }
    if (all) {
        if (large) for (i=ndims; i < SF_MAX_DIM; i++) {
	    sprintf(key,"n%d",i+1);
	    if (sf_histlargeint(in,key,&nlarge[i]) && nlarge[i]>0) ndims++;
	    else break; }
        else for (i=ndims; i < SF_MAX_DIM; i++) {
	    sprintf(key,"n%d",i+1);
	    if (sf_histint(in,key,&n[i]) && n[i]>0) ndims++;
	    else break; }
    }

    if (parform) printf( "ndims=" );

    printf( "%d", ndims );

    if (parform) printf( "\nn=" );
    else printf( ":" );
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

    if (parform) printf( "\nd=" );
    else printf( ":" );
    for (i=0; i < ndims; i++) {
        sprintf(key,"d%d",i+1);
        sf_histfloat(in,key,&fdummy);
	printf("%f",fdummy);
        if (i<ndims-1) printf(",");
    }

    if (parform) printf( "\no=" );
    else printf( ":" );
    for (i=0; i < ndims; i++) {
        sprintf(key,"o%d",i+1);
        sf_histfloat(in,key,&fdummy);
	printf("%f",fdummy);
        if (i<ndims-1) printf(",");
    }

    if (parform) printf( "\nlabel=" );
    else printf( ":" );
    for (i=0; i < ndims; i++) {
        sprintf(key,"label%d",i+1);
	printf("%s",sf_histstring(in,key)?sf_histstring(in,key):"none");
        if (i<ndims-1) printf(",");
    }

    if (parform) printf( "\nunit=" );
    else printf( ":" );
    for (i=0; i < ndims; i++) {
        sprintf(key,"unit%d",i+1);
	printf("%s",sf_histstring(in,key)?sf_histstring(in,key):"none");
        if (i<ndims-1) printf(",");
    }

    printf( "\n" );


    exit(0);
}
