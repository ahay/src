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
    float perc,fact, wp, *din=NULL, *dout=NULL;
    sf_complex *cin=NULL, *cout=NULL;
    sf_datatype type;
    sf_file in=NULL, out=NULL, other=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    n = sf_filesize(in);
    type = sf_gettype(in);

    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (!sf_getfloat("fact",&fact)) fact=0.5;
    /* factor for sharpening */

    sf_sharpen_init(n,perc,fact);

    if (SF_FLOAT == type) {

      din = sf_floatalloc(n);	
      dout = sf_floatalloc(n);

      sf_floatread(din,n,in);	

      if (NULL != sf_getstring("other")) {
	  other = sf_input("other");
	  if (SF_FLOAT != sf_gettype(other)) 
	      sf_error("Need float type in other");
	  sf_floatread(dout,n,other);
	  wp = sf_sharpen(dout);
      } else {
	  wp = sf_sharpen(din);
      }
      sf_warning("wp=%g",wp);

      sf_weight_lop(false,false,n,n,din,dout);
      sf_weight_lop(true,false,n,n,din,dout);

      sf_floatwrite(din,n,out);

    } else if (SF_COMPLEX == type) {

      cin = sf_complexalloc(n);	
      cout = sf_complexalloc(n);

      sf_complexread(cin,n,in);	

      if (NULL != sf_getstring("other")) {
	other = sf_input("other");
	if (SF_COMPLEX != sf_gettype(other)) 
	    sf_error("Need complex type in other");
	sf_complexread(cout,n,other);
	sf_csharpen(cout);
      } else {
	sf_csharpen(cin);
      }

      sf_cweight_lop(false,false,n,n,cin,cout);
      sf_cweight_lop(true,false,n,n,cin,cout);

      sf_complexwrite(cin,n,out);

    } else {
      sf_error("Need float or complex type");
    }

    exit(0);
}
