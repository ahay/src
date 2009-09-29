/* Image enhancement by histogram equalization. */
/*
  Copyright (C) 2006 University of Texas at Austin

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

#define NCOL 256

int main (int argc, char* argv[])
{
    int i3, n3, n1, n2, n12, i, col, accum;
    unsigned char* dat;
    int range[NCOL];
    sf_file in=NULL, out=NULL;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_UCHAR != sf_gettype(in)) sf_error("Need uchar type");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    dat = sf_ucharalloc(n12);

    for (i3=0; i3 < n3; i3++) {
	sf_ucharread(dat, n12, in);

	for (col=0; col < NCOL; col++) {
	    range[col] = 0;
	}

	for (i=0; i < n12; i++) {
	    col = dat[i];
	    range[col]++;
	}

	for (accum=0, col=0; col < NCOL; col++) {
	    accum += range[col];
	    range[col] = accum*(NCOL-1)/n12;
	}
  
	for (i=0; i < n12; i++) {
	    col = dat[i];
	    dat[i] = range[col];
	}

	sf_ucharwrite(dat, n12, out);
    }


    exit (0);
}
