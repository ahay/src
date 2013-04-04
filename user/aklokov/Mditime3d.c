/* 3D Hybrid Radon transform for diffraction imaging in the time dip-angle domain */
/*
  Copyright (C) 2013 University of Texas at Austin

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

#include "ditime3d.h"

int main (int argc, char* argv[]) {

    sf_file in, out, fileDweight, fileRefl, fileDiff;

	// data
    int tn; float to, td;
	int dipn;  float dipo, dipd;
    int sdipn; float sdipo, sdipd;    
	// reflection model
	int dip0n;  float dip0o, dip0d; // dips in x-direction
	int sdip0n; float sdip0o, sdip0d; // dips in y-direction
	// diffraction model
    int xin;  float xio, xid; // xi in x-direction
    int sxin; float sxio, sxid; // xi in y-direction

    int xn, yn, dagNum;
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

    if (!sf_histint   (in, "n1", &tn) ) sf_error ("No n1= in input");
    if (!sf_histfloat (in, "o1", &to) ) sf_error ("No o1= in input");
    if (!sf_histfloat (in, "d1", &td) ) sf_error ("No d1= in input");

    if (!sf_histint   (in, "n4", &xn) ) sf_error ("No n4= in input");
    if (!sf_histint   (in, "n5", &yn) ) sf_error ("No n5= in input");
	dagNum = xn * yn;
		

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

		if (!sf_histfloat (in, "o3", &sdipo) ) sf_error ("No o3= in input");
		if (!sf_histfloat (in, "d3", &sdipd) ) sf_error ("No d3= in input");
		if (!sf_histint   (in, "n3", &sdipn) ) sf_error ("No n3= in input");

		// reflection model
		if (!sf_getint ("dip0n", &dip0n) ) sf_error ("Need dip0n=");
		/* number of dip0 values (if adj=y) */
		if (!sf_getfloat ("dip0d", &dip0d)) sf_error("Need dip0d=");
		/* dip0 sampling (if adj=y) */
		if (!sf_getfloat("dip0o",&dip0o)) sf_error("Need dip0d0=");
		/* dip0 origin (if adj=y) */
		if (!sf_getint ("sdip0n", &sdip0n) ) sf_error ("Need sdip0n=");
		/* number of sdip0 values (if adj=y) */
		if (!sf_getfloat ("sdip0d", &sdip0d)) sf_error ("Need sdip0d=");
		/* sdip0 sampling (if adj=y) */
		if (!sf_getfloat("sdip0o", &sdip0o)) sf_error ("Need sdip0d0=");
		/* sdip0 origin (if adj=y) */

		// diffraction model
		if (!sf_getint ("xin", &xin) ) sf_error ("Need xin=");
		/* number of xi values (if adj=y) */
		if (!sf_getfloat ("xid", &xid)) sf_error("Need xid=");
		/* xi sampling (if adj=y) */
		if (!sf_getfloat("xio",&xio)) sf_error("Need xio=");
		/* xi origin (if adj=y) */
		if (!sf_getint ("sxin", &sxin) ) sf_error ("Need sxin=");
		/* number of xi values (if adj=y) */
		if (!sf_getfloat ("sxid", &sxid)) sf_error("Need sxid=");
		/* xi sampling (if adj=y) */
		if (!sf_getfloat("sxio", &sxio)) sf_error("Need sxio=");
		/* xi origin (if adj=y) */

		// diffraction-model file
		sf_putint    (fileDiff,  "n2", xin);
		sf_putfloat  (fileDiff,  "d2", xid);
		sf_putfloat  (fileDiff,  "o2", xio);
		sf_putstring (fileDiff, "label2", "xi");
		sf_putstring (fileDiff, "unit2", "");

		sf_putint    (fileDiff,  "n3", sxin);
		sf_putfloat  (fileDiff,  "d3", sxid);
		sf_putfloat  (fileDiff,  "o3", sxio);
		sf_putstring (fileDiff, "label3", "xi");
		sf_putstring (fileDiff, "unit3", "");

		// reflection-model file
		sf_putint    (fileRefl, "n2", dip0n);
		sf_putfloat  (fileRefl, "d2", dip0d);
		sf_putfloat  (fileRefl, "o2", dip0o);
		sf_putstring (fileRefl, "label2", "dip angle");
		sf_putstring (fileRefl, "unit2", "deg");

		sf_putint    (fileRefl, "n3", sdip0n);
		sf_putfloat  (fileRefl, "d3", sdip0d);
		sf_putfloat  (fileRefl, "o3", sdip0o);
		sf_putstring (fileRefl, "label3", "dip angle");
		sf_putstring (fileRefl, "unit3", "deg");

    } else { // model -> data

		fileDiff = in;
		fileRefl = sf_input ("reflMod");

		// reflection model
		if ( !sf_histint   (fileRefl, "n2", &dip0n) ) sf_error ("No n2= in reflection-model file");	
		if ( !sf_histfloat (fileRefl, "d2", &dip0d) ) sf_error ("No d2= in reflection-model file");	
		if ( !sf_histfloat (fileRefl, "o2", &dip0o) ) sf_error ("No o2= in reflection-model file");

		if ( !sf_histint   (fileRefl, "n3", &sdip0n) ) sf_error ("No n3= in reflection-model file");	
		if ( !sf_histfloat (fileRefl, "d3", &sdip0d) ) sf_error ("No d3= in reflection-model file");	
		if ( !sf_histfloat (fileRefl, "o3", &sdip0o) ) sf_error ("No o3= in reflection-model file");

		// diffraction model
		if ( !sf_histint   (fileDiff, "n2", &xin) ) sf_error ("No n2= in diffraction-model file");	
		if ( !sf_histfloat (fileDiff, "d2", &xid) ) sf_error ("No d2= in diffraction-model file");	
		if ( !sf_histfloat (fileDiff, "o2", &xio) ) sf_error ("No o2= in diffraction-model file");	

		if ( !sf_histint   (fileDiff, "n3", &sxin) ) sf_error ("No n3= in diffraction-model file");		
		if ( !sf_histfloat (fileDiff, "d3", &sxid) ) sf_error ("No d3= in diffraction-model file");		
		if ( !sf_histfloat (fileDiff, "o3", &sxio) ) sf_error ("No o3= in diffraction-model file");	

		// run parameters
		if (!sf_getint ("dipn", &dipn)) sf_error ("Need dipn=");
		/* number of dips in x-direction */
		if (!sf_getfloat ("dipo", &dipo)) sf_error ("Need dipo=");
		/* dip origin in x-direction */
		if (!sf_getfloat ("dipd", &dipd)) sf_error ("Need dipd=");
		/* dip sampling in x-direction */
		if (!sf_getint ("sdipn", &sdipn)) sf_error ("Need sdipn=");
		/* number of dips in y-direction */
		if (!sf_getfloat ("sdipo", &sdipo)) sf_error ("Need sdipo=");
		/* dip origin in y-direction */
		if (!sf_getfloat ("sdipd", &sdipd)) sf_error ("Need sdipd=");
		/* dip sampling in y-direction */

		// output file
		sf_putint    (out, "n2", dipn);
		sf_putfloat  (out, "o2", dipo);
		sf_putfloat  (out, "d2", dipd);
		sf_putstring (out, "label2", "x-dip angle");
		sf_putstring (out, "unit2", "deg");

		sf_putint    (out, "n3", sdipn);
		sf_putfloat  (out, "o3", sdipo);
		sf_putfloat  (out, "d3", sdipd);
		sf_putstring (out, "label3", "y-dip angle");
		sf_putstring (out, "unit3", "deg");
    }

	// data size
    dataSize = tn * dipn * sdipn;
    // model sizes
    diffSize = tn * xin * sxin; 
    reflSize = tn * dip0n * sdip0n;
    
	modelSize = -1;
    switch (invMod) {
		case 0: modelSize = diffSize; break;
		case 1: modelSize = reflSize; break;
		case 2: modelSize = diffSize + reflSize; break;
    }

    data    = sf_floatalloc (dataSize);
    model  = sf_floatalloc (modelSize);

    w = (0 == niter) ? NULL : sf_floatalloc (modelSize);
    p = (0 == niter) ? NULL : sf_floatalloc (modelSize);

    ditime3d_init (dipo, dipd, dipn, 
				   sdipo, sdipd, sdipn, 
				   xio, xid, xin, 
				   sxio, sxid, sxin, 
				   dip0o, dip0d, dip0n, 
				   sdip0o, sdip0d, sdip0n, 
				   to, td, tn, 
				   isAA, invMod);

	// main loop	
	
    for (ix = 0; ix < dagNum; ++ix) { 
		if (verb) sf_warning ("i=%d of %d", ix + 1, dagNum);

		// read data	
		if (adj) { // data -> model
		    sf_floatread (data, dataSize, in);
		    dweight = sf_floatalloc (dataSize);
		    if (fileDweight) {
				sf_floatread (dweight, dataSize, fileDweight);
			} else {
				for (id = 0; id < dataSize; ++id) dweight[id] = 1.f;
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
		    ditime3d_lop (adj, false, modelSize, dataSize, model, data);
		} else {
		    // initital model weights
		    for (im = 0; im < modelSize; ++im) {
				w[im] = 1.f;
		    }

		    for (iter = 0; iter < niter; ++iter) {
				sf_solver_prec (ditime3d_lop, sf_cgstep, sf_copy_lop,
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
		} else { // model -> data
		    sf_floatwrite (data, dataSize, out);
		} 
    }

	// finish
	ditime3d_close ();

	if (w) free (w);
	if (p) free (p);
	if (dweight) free (dweight);    

    free (data);
    free (model);

    sf_fileclose (in);
    sf_fileclose (out);
    sf_fileclose (fileRefl);
    if (fileDweight) sf_fileclose (fileDweight);

    return 0;
}
