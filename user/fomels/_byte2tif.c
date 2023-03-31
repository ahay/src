/* Convert byte RSF to a JPEG image. */
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

#include "_jpeg.h"

int main(int argc, char* argv[])
{
    bool color;
    int n1, n2, nc;
    unsigned char *grey=NULL;
    sf_file in=NULL;

    sf_init(argc,argv);
    in = sf_input("in");

    if (SF_UCHAR != sf_gettype(in)) sf_error("Need unsigned char in input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_getbool("color",&color)) color=(bool)(3==n1);

    if (color) {
	nc = n1;
	if (!sf_histint(in,"n2",&n1)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n2)) sf_error("No n3= in input");
    } else {
	nc = 1;
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    }

    grey = sf_ucharalloc (n1*n2*nc);
    sf_ucharread(grey,n1*n2*nc,in);    
    write_JPEG_file (grey, n2, n1, nc);

    exit(0);
}
