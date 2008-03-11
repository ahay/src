/* Sharpening operator */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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
    int n;
    float perc,a, *din, *dout;
    sf_file in, out, other;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    n = sf_filesize(in);
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
     
    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (!sf_getfloat("perc",&a)) a=0.1;
    /* sharpening cutoff */

    sf_sharpen_init(n,perc);

    din = sf_floatalloc(n);
    dout = sf_floatalloc(n);

    sf_floatread(din,n,in);

    if (NULL != sf_getstring("other")) {
	other = sf_input("other");
	sf_floatread(dout,n,other);
	sf_sharpen(dout);
    } else {
	sf_sharpen(din);
    }

    sf_weight_lop(false,false,n,n,din,dout);
    sf_weight_lop(true,false,n,n,din,dout);

    sf_floatwrite(din,n,out);

    exit(0);
}
