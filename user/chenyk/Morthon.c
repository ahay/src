/* Non-stationary orthogonalization */
/*
  Copyright (C) 2020 Zhejiang University & University of Texas at Austin
   
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
#include "ntrianglen.h"
#include "divnn.h"

int main(int argc, char* argv[])
{
    bool verb;
    int   *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */
    int i, id, dim, n[SF_MAX_DIM], dim1, nd, niter, n1, n2, i1, b;
    float *noi, *sig, *rat, *noi2, *sig2, eps, remove;
    char key[8];
    int box[SF_MAX_DIM];
	/*box means maximum (edge) padding for triangle smoothing, box[i]=max(rect[i])*/
    sf_file fnoi, fsig, fnoi2, fsig2;
	sf_file rect[SF_MAX_DIM], shift[SF_MAX_DIM]; 
    
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
	nd *= n[i];
    }
    
	/*Calculate dim1*/
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    rect[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
	    dim1 = i;
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rect[i] = NULL;
	    shift[i] = NULL;
	}
    }
    
	/*Calculate n1*/
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }
    
    noi = sf_floatalloc(nd);
    sig = sf_floatalloc(nd);
    rat = sf_floatalloc(nd);    

    noi2 = sf_floatalloc(nd);
    sig2 = sf_floatalloc(nd);

	/*reading the non-stationary smoothing radii*/
    for (i=0; i <= dim1; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_floatalloc (n1);
	    sft[i] = sf_intalloc (n1);

	    sf_floatread(rct[i],n1,rect[i]);
	    sf_fileclose(rect[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[i],n1,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    sft[i][i1] = 0;
		}
	    }

		
	    for (i1=0; i1 < n1; i1++) {
		b = ceilf(rct[i][i1])+SF_ABS(sft[i][i1]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

	sf_warning("dim=%d, dim1=%d, nd=%d, n1=%d",dim, dim1, nd,n1);
	sf_warning("n[0]=%d, n[1]=%d, n[2]=%d",n[0],n[1],n[2]);	
	sf_warning("box[0]=%d, box[1]=%d, box[2]=%d",box[0],box[1],box[2]);	
	/*box means maximum (edge) padding for triangle smoothing, box[i]=max(rect[i])*/
	sf_warning("nd=%d",nd);
	
    divnn_init(dim, nd, box, n, rct, sft, niter, verb);
    
    sf_floatread(noi,nd,fnoi);
    sf_floatread(sig,nd,fsig);

    for (id=0; id < nd; id++) {
	noi2[id] = noi[id];
	sig2[id] = sig[id];
    }

    divnne (noi, sig, rat, eps);
    
    for (id=0; id < nd; id++) {
	remove = rat[id]*sig2[id];
	noi2[id] -= remove;
	sig2[id] += remove;
    }

    sf_floatwrite(noi2,nd,fnoi2);
    sf_floatwrite(sig2,nd,fsig2);

    exit(0);
}
