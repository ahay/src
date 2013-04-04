/* 2D Hybrid Radon transform for diffraction imaging in the time dip-angle domain */
/*
  Copyright (C) 2012 University of Texas at Austin

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

#include "ditime2d.h"

int main (int argc, char* argv[]) {

    sf_file in, out, fileDweight, fileRefl, fileDiff;
    int dipn;  float dipo, dipd;
    int dip0n; float dip0o, dip0d;
    int tn; float to, td;
    int xin;   float xio, xid;
    int xn;
    float *data = NULL, *model = NULL, *dweight = NULL;
    int im, ix;
    int invMod;

    int   niter, liter, iter;
    bool  adj, verb, isAA;
    float eps;
    float *w, *p;

    int dataSize, diffSize, reflSize, modelSize, id;

    sf_init (argc, argv);

    in  = sf_input  ("in");
    out = sf_output ("out");
	
    fileDweight = NULL;
    if ( NULL != sf_getstring ("dweight") ) {
	/* input file containing data weights */ 
	fileDweight = sf_input ("dweight");
    }

    if (!sf_histint   (in, "n1", &tn) ) sf_error("No n1= in input");
    if (!sf_histfloat (in, "o1", &to) ) sf_error("No o1= in input");
    if (!sf_histfloat (in, "d1", &td) ) sf_error("No d1= in input");

    xn = sf_leftsize (in, 2);

    if (!sf_getbool ("verb", &verb)) verb=false; /* verbosity flag */
    if (!sf_getbool ( "adj", &adj )) adj=false;  /* adjoint flag */
    if (!sf_getbool ("isAA", &isAA)) isAA = false;
    /* if y, apply anti-aliasing */

    if (!sf_getint ("liter",&liter)) liter=100; /* number of linear iterations (for inversion) */
    if (!sf_getint ("niter",&niter)) niter=0; /* number of nonlinear iterations (for inversion) */
    if (!sf_getfloat ("eps",&eps)) eps=0.;    /* regularization parameter */
    if (!sf_getbool("verb",&verb)) verb = false; /* verbosity flag */
    if (!sf_getint ("invMod", &invMod)) invMod=2; /* number of nonlinear iterations (for inversion) */

    if (adj) { // data -> model
		fileDiff = out;
		fileRefl = sf_output ("reflMod");

		if (!sf_histfloat (in, "o2", &dipo) ) sf_error ("No o2= in input");
		if (!sf_histfloat (in, "d2", &dipd) ) sf_error ("No d2= in input");
		if (!sf_histint   (in, "n2", &dipn) ) sf_error ("No n2= in input");

		// REFLECTION MODEL
		if (!sf_getint ("dip0n", &dip0n) ) sf_error ("Need dip0n=");
		/* number of dip0 values (if adj=y) */
		if (!sf_getfloat ("dip0d", &dip0d)) sf_error("Need dip0d=");
		/* dip0 sampling (if adj=y) */
		if (!sf_getfloat("dip0o",&dip0o)) sf_error("Need dip0d0=");
		/* dip0 origin (if adj=y) */

		// DIFFRACTION MODEL
		if (!sf_getint ("xin", &xin) ) sf_error ("Need xin=");
		/* number of xi values (if adj=y) */
		if (!sf_getfloat ("xid", &xid)) sf_error("Need xid=");
		/* xi sampling (if adj=y) */
		if (!sf_getfloat("xio",&xio)) sf_error("Need xio=");
		/* xi origin (if adj=y) */

		sf_putint   (fileDiff,  "n2", xin);
		sf_putfloat (fileDiff,  "d2", xid);
		sf_putfloat (fileDiff,  "o2", xio);
		sf_putstring (fileDiff, "label2", "xi");
		sf_putstring (fileDiff, "unit2", "");

		sf_putint   (fileRefl,  "n2", dip0n);
		sf_putfloat (fileRefl,  "d2", dip0d);
		sf_putfloat (fileRefl,  "o2", dip0o);
		sf_putstring (fileRefl, "label2", "dip angle");
		sf_putstring (fileRefl, "unit2", "deg");

    } else { // model -> data

		fileDiff = in;
		fileRefl = sf_input ("reflMod");

		if (!sf_histint   (fileDiff, "n2", &xin)) sf_error("No n2= in input");	
		if (!sf_histfloat (fileDiff, "d2", &xid)) sf_error("No d2= in input");	
		if (!sf_histfloat (fileDiff, "o2", &xio)) sf_error("No o2= in input");
	
		if (!sf_histint   (fileRefl, "n2", &dip0n)) sf_error("No n2= in input");	
		if (!sf_histfloat (fileRefl, "d2", &dip0d)) sf_error("No d2= in input");	
		if (!sf_histfloat (fileRefl, "o2", &dip0o)) sf_error("No o2= in input");

		if (!sf_getint ("dipn", &dipn)) sf_error ("Need dipn=");
		/* number of offsets */
		if (!sf_getfloat ("dipo", &dipo)) sf_error("Need dipo=");
		/* offset origin */
		if (!sf_getfloat ("dipd", &dipd)) sf_error("Need dipd=");
		/* offset sampling */

		sf_putint    (out, "n2", dipn);
		sf_putfloat  (out, "o2", dipo);
		sf_putfloat  (out, "d2", dipd);
		sf_putstring (out, "label2", "dip angle");
		sf_putstring (out, "unit2", "deg");
    }

	// data size
    dataSize = tn * dipn;
    // model sizes
    diffSize = tn * xin;
    reflSize = tn * dip0n;

    modelSize = -1;
    switch (invMod) {
		case 0: modelSize = diffSize; break;
		case 1: modelSize = reflSize; break;
		case 2: modelSize = diffSize + reflSize; break;
    }

    data   = sf_floatalloc (dataSize);
    model  = sf_floatalloc (modelSize);

    w = (0 == niter) ? NULL : sf_floatalloc (modelSize);
    p = (0 == niter) ? NULL : sf_floatalloc (modelSize);

    ditime2d_init (dipo, dipd, dipn, 
				   xio, xid, xin, 
				   dip0o, dip0d, dip0n, 
				   to, td, tn, 
				   isAA, invMod);

    for (ix = 0; ix < xn; ++ix) { 
		if (verb) sf_warning ("i=%d of %d", ix + 1, xn);

		// read data	
		if (adj) { // data -> model
		    sf_floatread (data, dataSize, in);
		    dweight = sf_floatalloc (dataSize);
		    if (fileDweight) {
				sf_floatread (dweight, dataSize, fileDweight);
			} else {
				for (id = 0; id < dataSize; ++id)
				    dweight[id] = 1.f;
			}
		} else { // model -> data
			if (0 == invMod) { // difffraction
			    sf_floatread (model, diffSize, fileDiff);
			} else if (1 == invMod) { // reflection
				sf_floatread (model, reflSize, fileRefl);
			} else { // diffraction + reflection
			    sf_floatread (model, diffSize, fileDiff);
			    sf_floatread (model + diffSize, reflSize, fileRefl);
			}
		}
		
		// perform transform
		if (!adj || 0 == niter) {
		    ditime2d_lop (adj, false, modelSize, dataSize, model, data);
		} else {
		    // initital model weights
		    for (im = 0; im < modelSize; ++im) {
				w[im] = 1.f;
		    }

		    for (iter = 0; iter < niter; ++iter) {
				sf_solver_prec (ditime2d_lop, sf_cgstep, sf_copy_lop,
								modelSize, modelSize, dataSize, model, data, liter, eps,
								"verb", verb, "mwt", w, "xp", p, "wt", dweight, "end");
				sf_cgstep_close ();

				for (im = 0; im < modelSize; ++im) {
				    w[im] = fabsf (p[im]); /* weight for sparsity */
				}
		    }
		}
		// write result
		if (adj) { // data -> model
			if (0 == invMod) { // difffraction
			    sf_floatwrite (model, diffSize, fileDiff);
			} else if (1 == invMod) { // reflection
				sf_floatwrite (model, reflSize, fileRefl);
			} else { // diffraction + reflection
			    sf_floatwrite (model, diffSize, fileDiff);
			    sf_floatwrite (model + diffSize, reflSize, fileRefl);
			}
		} else {
		    sf_floatwrite (data, dataSize, out);
		} 
    }

	// finish

	ditime2d_close ();

    sf_fileclose (in);
    sf_fileclose (out);
    sf_fileclose (fileRefl);
    if (fileDweight) sf_fileclose (fileDweight);

    free (data);
    free (model);

    if (dweight) free (dweight);
	if (w) free (w);
	if (p) free (p);

    return 0;
}
