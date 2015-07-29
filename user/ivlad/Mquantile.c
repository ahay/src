/* Computes what clip value corresponds to a given pclip.
Loads the entire dataset in core. Use it to find a clip= parameter for sfclip, given a wanted pclip=*/
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
    int n; /* Total number of samples in dataset */
    float pclip;
    int nclip;
    float q; /* Quantile */
    float *d=NULL; /* Array to hold THE ENTIRE DATASET */
    sf_file in=NULL; /* Input file */
    int mem; /* For avoiding int to off_t typecast warning */
    off_t memsize;

    sf_init(argc,argv);

    if( !sf_getfloat("pclip", &pclip))
        sf_error( "Need pclip=" );
        /* Percentage clip, between 0 and 100 */

    /* Give the user high-resolution error messages */
    if( pclip<0  ) sf_error( "pclip must be >0"   );
    if( pclip>100) sf_error( "pclip must be <100" );

    in = sf_input("in");

    if (SF_FLOAT != sf_gettype(in))
        sf_error( "Need float input" );

    if( !sf_getint("memsize",&mem))
        mem=sf_memsize();
    /* Max amount of RAM (in Mb) to be used */
    memsize = mem * (1<<20); /* convert Mb to bytes */

    n = sf_leftsize(in,0);

    if( memsize < n*sizeof(float) )
        sf_error( "Not enough memory" );

    d = sf_floatalloc(n);

    sf_floatread( d, n, in );

    /* Shoplifted from plot/lib/gainpar.c: */
    nclip = SF_MAX(SF_MIN(n*pclip/100. + .5,n-1),0);

    q = sf_quantile( nclip, n, d );

    printf( "%f\n", q );


    exit(0);
}
