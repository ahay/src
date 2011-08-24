/* 1-D convolution with complex numbers. */
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

#include "cicai1.h"

int main(int argc, char* argv[])
{
    int nx, nf, i2, n2, lag;
    bool single;
    sf_complex *xx, *yy, *ff;
    sf_file in, out, filt;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    filt = sf_input("filt");

    if (SF_COMPLEX != sf_gettype(in) ||
	SF_COMPLEX != sf_gettype(filt)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(filt,"n1",&nf)) sf_error("No n1= in filtin");

    xx = sf_complexalloc(nx);
    ff = sf_complexalloc(nf);

    if (!sf_getbool("single",&single)) single=true;
    /* single channel or multichannel */

    if (!sf_getint("lag",&lag)) lag=1;
    /* lag for internal convolution */

    yy = sf_complexalloc(nx);

    if (!single) sf_complexread (ff,nf,filt);

    for (i2=0; i2 < n2; i2++) {
        if (single) sf_complexread (ff,nf,filt);
        sf_complexread (xx,nx,in);

	cicai1_init(nf,ff,lag);
	cicai1_lop (false, false,nx,nx,xx,yy);
	
	sf_complexwrite (yy,nx,out);
    }
    
    exit(0);
}

/* 	$Id: Mconv.c 7107 2011-04-10 02:04:14Z ivlad $	 */
