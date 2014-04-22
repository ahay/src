/* Orthogonolize signal and noise. */
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
    int i, id, dim, n[SF_MAX_DIM], nd, rect[SF_MAX_DIM], niter;
    float *noi, *sig, *rat, *noi2, *sig2, eps, remove;
    char key[6];
    sf_file fnoi, fsig, fnoi2, fsig2;

    sf_init(argc,argv);
    fnoi = sf_input("in");
    fsig = sf_input("sig"); /* input signal */

    if (SF_FLOAT != sf_gettype(fnoi) ||
	SF_FLOAT != sf_gettype(fsig)) sf_error("Need float input");

    fnoi2 = sf_output("out");
    fsig2 = sf_output("sig2"); /* output signal */

    dim = sf_filedims (fnoi,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	nd *= n[i];
    }

    noi = sf_floatalloc(nd);
    sig = sf_floatalloc(nd);
    rat = sf_floatalloc(nd);    

    noi2 = sf_floatalloc(nd);
    sig2 = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    sf_divn_init(dim, nd, n, rect, niter, verb);

    sf_floatread(noi,nd,fnoi);
    sf_floatread(sig,nd,fsig);

    for (id=0; id < nd; id++) {
	noi2[id] = noi[id];
	sig2[id] = sig[id];
    }

    sf_divne (noi, sig, rat, eps);

    for (id=0; id < nd; id++) {
	remove = rat[id]*sig2[id];
	noi2[id] -= remove;
	sig2[id] += remove;
    }

    sf_floatwrite(noi2,nd,fnoi2);
    sf_floatwrite(sig2,nd,fsig2);

    exit(0);
}
