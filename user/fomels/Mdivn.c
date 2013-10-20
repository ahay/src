/* Smooth division. */
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

int main(int argc, char* argv[])
{
    bool verb;
    int i, dim, n[SF_MAX_DIM], nd, rect[SF_MAX_DIM], niter;
    float *num, *den, *rat, eps;
    char key[6];
    sf_file fnum, fden, frat;

    sf_init(argc,argv);
    fnum = sf_input("in");
    fden = sf_input("den");
    frat = sf_output("out");

    if (SF_FLOAT != sf_gettype(fnum) ||
	SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");

    dim = sf_filedims (fnum,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	nd *= n[i];
    }

    num = sf_floatalloc(nd);
    den = sf_floatalloc(nd);
    rat = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    sf_divn_init(dim, nd, n, rect, niter, verb);

    sf_floatread(num,nd,fnum);
    sf_floatread(den,nd,fden);

    sf_divne (num, den, rat, eps);

    sf_floatwrite(rat,nd,frat);

    exit(0);
}
