/* Computes Ni+1 x Ni+2 x ...
Wrapper for sf_leftsize*/
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
    int n;
    int i;
    sf_file in=NULL;

    sf_init(argc,argv);

    if( !sf_getint("i", &i)) i=0;
    /* What size to start counting from. i=0 gets total number of elements */

    if( i < 0  ) sf_error( "i must be >= 0" );

    in = sf_input("in");

    if (SF_FLOAT != sf_gettype(in))
        sf_error( "Need float input" );

    n = sf_leftsize(in,i);

    printf( "%d\n", n );


    exit(0);
}
