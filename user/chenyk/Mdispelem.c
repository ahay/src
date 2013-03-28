/* Display element of rsf files. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

int main (int argc, char *argv[])
{
/*
    bool verb;
*/
    int n1, n2, n12, i1,i2;
    float *trace;
    sf_file in;

    sf_init(argc,argv);
    in  = sf_input("in");
    /*out = sf_output("out");*/

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12=n1*n2;


    if (!sf_getint("i1",&i1)) sf_error("Need i1=");
    /* get the index of first axis */

    if (!sf_getint("i2",&i2)) sf_error("Need i2=");
    /* get the index of second axis */
  
    trace = sf_floatalloc(n12);
 
    sf_floatread(trace,n12,in);
    sf_warning("The element is %f", trace[(i1-1)+(i2-1)*n1]);
    
    exit (0);
}

