/* Smooth division for complex data. */
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

#include "cdivn.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, id, dim, n[SF_MAX_DIM], nd, rect[SF_MAX_DIM], niter;
    sf_complex *num, *den, *rat;
    float a, norm;
    char key[6];
    sf_file fnum, fden, frat;

    sf_init(argc,argv);
    fnum = sf_input("in");
    fden = sf_input("den");
    frat = sf_output("out");

    if (SF_COMPLEX != sf_gettype(fnum) ||
	SF_COMPLEX != sf_gettype(fden)) sf_error("Need complex input");

    dim = sf_filedims (fnum,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	nd *= n[i];
    }

    num = sf_complexalloc(nd);
    den = sf_complexalloc(nd);
    rat = sf_complexalloc(nd);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    cdivn_init(dim, nd, n, rect, niter, verb);

    sf_complexread(num,nd,fnum);
    sf_complexread(den,nd,fden);

    norm = 0.;
    for (id=0; id < nd; id++) {
	a = cabsf(den[id]);
	norm += a*a;
    }
    norm = sqrtf(nd/norm);

    for (id=0; id < nd; id++) {
#ifdef SF_HAS_COMPLEX_H
	num[id] *= norm;
	den[id] *= norm;
#else
	num[id] = sf_crmul(num[id],norm);
	den[id] = sf_crmul(den[id],norm);
#endif
    }

    cdivn (num, den, rat);

    sf_complexwrite(rat,nd,frat);

    exit(0);
}
