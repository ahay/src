/* Focusing indicator. 

Takes: rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
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
#include <math.h>

#include <rsf.h>

#include "divn.h"
#include "hilbert.h"

int main (int argc, char* argv[])
{
    int i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *num, *den, *rat, mean;
    char key[6];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	n12 *= n[i];
    }

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    rat = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    sf_floatread(den,n12,in);

    mean=0.;
    for (i=0; i < n12; i++) {
	mean += den[i]*den[i];
    }
    mean = sqrtf(mean/n12);

    for (i=0; i < n12; i++) {
	num[i] = den[i]*den[i]*den[i]*mean*n12;
	den[i] /= mean;
    }
    
    divn_init(dim, n12, n, rect, niter);
    divn (num, den, rat);

    sf_floatwrite(rat,n12,out);

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
