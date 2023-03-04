/* Smooth division. 

December 2015 program of the month:
http://ahay.org/blog/2015/12/22/program-of-the-month-sfdivn/
*/
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
    int i, dim, dim1, n[SF_MAX_DIM], rect[SF_MAX_DIM], niter, n1, i2, n2;
    float *num, *den, *rat, eps;
    char key[7];
    sf_file fnum, fden, frat;

    sf_init(argc,argv);
    fnum = sf_input("in");
    fden = sf_input("den");
    frat = sf_output("out");

    if (SF_FLOAT != sf_gettype(fnum) ||
	SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");

    dim = sf_filedims (fnum,n);

    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,7,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i+1;
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
        if (i < dim1) {
            n1 *= n[i];
        } else {
            n2 *= n[i];
        }
    }

    num = sf_floatalloc(n1);
    den = sf_floatalloc(n1);
    rat = sf_floatalloc(n1);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    sf_divn_init(dim1, n1, n, rect, niter, verb);

    for (i2=0; i2 < n2; i2++) {
	sf_warning(verb? 
		   "record %d of %d":
		   "record %d of %d;",i2+1,n2);

	sf_floatread(num,n1,fnum);
	sf_floatread(den,n1,fden);

	sf_divne (num, den, rat, eps);
	sf_floatwrite(rat,n1,frat);
    }
    if (!verb) sf_warning(".");

    exit(0);
}
