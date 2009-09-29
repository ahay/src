/* Local correlation with the envelope. */
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

#include <string.h>

#include <rsf.h> 

#include "envcorr.h"

int main(int argc, char* argv[])
{ 
    int dim, m[SF_MAX_DIM], nd, i, niter, rect[SF_MAX_DIM];
    float **inp, *rat;
    char key[6];
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    dim = sf_filedims(in,m);

    nd=1;
    for (i=0; i < dim; i++) {
	nd *= m[i];
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    
    inp = sf_floatalloc2 (m[0],nd/m[0]);
    rat = sf_floatalloc (nd);

    envcorr_init(dim,m,rect,niter);
    
    sf_floatread(inp[0],nd,in);

    envcorr(inp,rat);

    sf_floatwrite(rat,nd,out);

    exit (0);
}

/* 	$Id: Mwarpscan.c 744 2004-08-17 18:46:07Z fomels $	 */
