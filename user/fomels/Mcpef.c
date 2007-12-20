/* 1-D prediction-error filter estimation from complex data */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "cburg.h"

int main(int argc, char* argv[])
{
    int n1, i2, n2, nf;
    sf_complex *trace, *a;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("nf",&nf)) sf_error("Need nf=");
    /* filter length */
    sf_putint(out,"n1",nf);

    cburg_init(n1,1,nf);

    trace = sf_complexalloc(n1);
    a = sf_complexalloc(nf);

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(trace,n1,in);

	cburg_apply(trace,a);

	sf_complexwrite(a,nf,out);
    }

    exit(0);
}
