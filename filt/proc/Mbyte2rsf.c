/* Convert raw byte images to RSF.

Takes: < raw.img
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

#include <stdio.h>

#include <rsf.h> 

int main(int argc, char* argv[])
{ 
    int n1, n2, x, y;
    unsigned char *line;
    float* array;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");

    if(!sf_getint("n1",&n1)) sf_error("Need n1="); 
    /* vertical dimension */
    if(!sf_getint("n2",&n2)) sf_error("Need n2="); 
    /* horizontal dimension */

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",1.);
    sf_putfloat(out,"d2",1.);
    sf_setformat(out,"native_float");

    array = sf_floatalloc(n1);
    line = (unsigned char *) sf_alloc(n1,sizeof(unsigned char));

    for (y = 0; y < n2; y++) {
	if (n1 != fread(line, sizeof(unsigned char), n1, stdin))
	    sf_error("trouble reading input data");
	for (x = 0; x < n1; x++) {
	    array[x] = (float) line[x];
	}
	sf_floatwrite(array,n1,out); 
    }

    exit (0);
}

/* 	$Id$	 */
