/* Argument of complex data calculated by atan2. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include <math.h>
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1,n2,i,i2;
    float x,y,*a;
    sf_complex *z;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file */
    sf_settype(in,SF_COMPLEX);
    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    /* output real */
    sf_settype(out,SF_FLOAT);

    /* memory allocations */
    z = sf_complexalloc(n1);
    a = sf_floatalloc(n1);

    for (i2 = 0; i2 < n2; i2++) {

	sf_complexread(z,n1,in);

	for (i = 0; i < n1; i++) {

	    x = crealf(z[i]);
	    y = cimagf(z[i]);

	    a[i] = atan2(y,x)*180.0/SF_PI;

	}

	sf_floatwrite(a,n1,out);
    }

    exit (0);
}


