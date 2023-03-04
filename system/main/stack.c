/* Stack a dataset over one of the dimensions.

   This operation is adjoint to sfspray.
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2021 Colorado School of Mines

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
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int j, ni, ntmp;
  
    int axis, ndims[SF_MAX_DIM],dim,lim;
    sf_file in, out;

    char key[7], *prog;
    bool norm, rms, min=false, max=false, prod=false, all=false;

    float          tR=0.0;
    float      * dinR;
    float      * douR;
    double     * stkR;

    sf_complex     tC=0.0;
    sf_complex * dinC;
    sf_complex * douC;
    sf_double_complex * stkC;

    float      * scale = NULL;
    int        * fold  = NULL;

    bool isreal; /* true for float; false for complex */

    off_t nBELOW, nAXIS, nABOVE;
    off_t iBELOW, iAXIS, iABOVE;

    sf_init (argc, argv);
    in  = sf_input ( "in");
    out = sf_output("out");

    /* ------------------------------------------------------------ */
    if      ( sf_gettype(in) == SF_FLOAT ) {
	isreal = true;
    } else if ( sf_gettype(in) == SF_COMPLEX ) {
	isreal = false;
    } else {
	isreal = false;
	sf_error("Incorrect data type in input");
    }
    /* ------------------------------------------------------------ */
  
    dim = sf_filedims(in,ndims);

    if (!sf_getint("axis",&axis)) axis = 2;
    /* which axis to stack. If axis=0, stack over all dimensions */

    if( axis == 0) {
	axis = dim + 1;
	all  = true;
	lim  = axis;
    } else {
	lim  = axis - 1;
    }

    /*
      sf_warning("dim  = %d",dim);
      sf_warning("axis = %d",axis);
      sf_warning("lim  = %d",lim);
    */

    nBELOW = 1;
    for ( j = 0; j < lim; j++) {
	snprintf(key,3,"n%d",(j+1)%10u);
	if (!sf_histint(in,key,&ni)) break;
	nBELOW *= ni;
    }

    /* ------------------------------------------------------------ */

    if(all) {
	nAXIS  = nBELOW;
	nABOVE = 1;
	nBELOW = 1;

	for (j = 0; j < lim; j++) {
	    snprintf(key,3,"n%d",(j+1)%10u);
	    sf_putint(out,key,1);
	}

    } else {

	snprintf(key,3,"n%d",axis%10u);
	if (!sf_histint(in,key,&ntmp)) sf_error("No %s= in input",key);
	nAXIS = ntmp;

	nABOVE = sf_unshiftdim(in,out,axis);
    }

    /*
      sf_warning("-------------------");
      sf_warning("nBELOW = %d",nBELOW ); // samples below the stack axis
      sf_warning("nAXIS  = %d",nAXIS);   // samples of the axis to stack
      sf_warning("nABOVE = %d",nABOVE);  // samples above the stack axis
      sf_warning("-------------------");
    */

    /* ------------------------------------------------------------ */

    if (nAXIS < 100) {
  	scale = sf_floatalloc(nAXIS);
	if (!sf_getfloats("scale",scale,nAXIS)) { // scale values on stack axis
	    /* optionally scale before stacking */
	    free(scale);
	    scale = NULL;
	}
    }

    /* ------------------------------------------------------------ */
    prog = sf_getprog();
    if        (NULL != strstr(prog,"max")) {
  	max = true;
    } else if (NULL != strstr(prog,"min")) {
  	min = true;
    } else if (NULL != strstr(prog,"prod")) {
  	prod = true;
    } else {

  	if (       !sf_getbool("rms" ,&rms))  rms  = false;
	/* If y, compute the root-mean-square instead of stack. */
	if (rms || !sf_getbool("norm",&norm)) norm = true;
	/* If y, normalize by fold. */
	if (       !sf_getbool("min" ,&min))  min  = false;
	/* If y, find minimum instead of stack. Ignores rms and norm. */
	if (       !sf_getbool("max" ,&max))  max  = false;
	/* If y, find maximum instead of stack. Ignores rms and norm. */
  	if (       !sf_getbool("prod",&prod)) prod = false;
  	/* If y, find product instead of stack. Ignores rms and norm. */

  	if ((min && max) || (min && prod) || (max && prod))
	    sf_error("Cannot have min=y max=y prod=y");

    }

    if (min || max || prod) { rms = false; norm = false; }
    /* ------------------------------------------------------------ */

    if(norm) fold = sf_intalloc(nBELOW);

    if(isreal) {
	dinR = sf_floatalloc(nBELOW);
	douR = sf_floatalloc(nBELOW);
	dinC = NULL;
	douC = NULL;
	stkR = (min || max) ? NULL :            (double*) sf_alloc(nBELOW,sizeof(double));
	stkC = NULL;
    } else {
	dinR = NULL;
	douR = NULL;
	dinC = sf_complexalloc(nBELOW);
	douC = sf_complexalloc(nBELOW);
	stkR = NULL;
	stkC = (min || max) ? NULL : (sf_double_complex*) sf_alloc(nBELOW,sizeof(sf_double_complex));
    }

    /* ------------------------------------------------------------ */

    for (iABOVE = 0; iABOVE < nABOVE; iABOVE++) {

	if(min || max) { // prepare to search for min or max
	    if(isreal) {
		if(min) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) douR[iBELOW] = +FLT_MAX;
		if(max) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) douR[iBELOW] = -FLT_MAX;
	    } else {
#ifdef SF_HAS_COMPLEX_H
		if(min) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) {
			douC[iBELOW] = +FLT_MAX+FLT_MAX*I;
		    }
		if(max) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) {
			douC[iBELOW] = -FLT_MAX-FLT_MAX*I;
		    }
#else
		if(min) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) {
			douC[iBELOW].r = +FLT_MAX;
			douC[iBELOW].i = +FLT_MAX;
		    }
		if(max) for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) {
			douC[iBELOW].r = -FLT_MAX;
			douC[iBELOW].i = -FLT_MAX;
		    }
#endif
	    }

	} else {

	    if(isreal) {
		for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) { // initialize prod w/ 1 or stack w/ 0
		    stkR[iBELOW] = prod ? 1.0 : 0.0;
		}
	    } else {
		for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) { // initialize prod w/ 1 or stack w/ 0
		    stkC[iBELOW] = prod ? 1.0 : 0.0;
		}
	    }

	}

	if(norm) memset (fold , 0 , nBELOW * sizeof(int));

	for (iAXIS = 0; iAXIS < nAXIS; iAXIS++) {
	    if(isreal) sf_floatread  (dinR, nBELOW, in);
	    else       sf_complexread(dinC, nBELOW, in);

	    for( iBELOW = 0; iBELOW < nBELOW; iBELOW++) {

		if(isreal) tR = dinR[iBELOW];
		else       tC = dinC[iBELOW];

		if (NULL != scale) {
		    if(isreal) tR *= scale[iAXIS];
		    else       tC *= scale[iAXIS];
		}

		if (min || max) {
          
		    if(isreal) {
			if(min) douR[iBELOW] = tR < douR[iBELOW] ? tR : douR[iBELOW];
			if(max) douR[iBELOW] = tR > douR[iBELOW] ? tR : douR[iBELOW];
		    } else {
			if(min) douC[iBELOW] = cabsf(tC) < cabsf(douC[iBELOW]) ? tC : douC[iBELOW];
			if(max) douC[iBELOW] = cabsf(tC) > cabsf(douC[iBELOW]) ? tC : douC[iBELOW];
		    }

		} else if (prod) {

		    if(isreal) stkR[iBELOW] *= tR;
		    else       stkC[iBELOW] *= tC;
            
		} else {

		    if(isreal) stkR[iBELOW] += rms ? (double)       tR*tR : tR;
		    else       stkC[iBELOW] += rms ? (double) conj(tC)*tC : tC;
            
		}

		if(isreal) { if (norm && (0.0 !=      tR))  fold[iBELOW]++; }
		else       { if (norm && (0.0 != cabs(tC))) fold[iBELOW]++; }
      
	    }
	}

	if ( !min && !max ) {

	    for (iBELOW = 0; iBELOW < nBELOW; iBELOW++) {

		if (norm) {
		    if (fold[iBELOW] > 0) {

			if(isreal) stkR[iBELOW] /= fold[iBELOW];
			else       stkC[iBELOW] /= fold[iBELOW];

			if(isreal) { if (rms) stkR[iBELOW] =  sqrt(stkR[iBELOW]); }
			else       { if (rms) stkC[iBELOW] = csqrt(stkC[iBELOW]); }
		    }
		}

		if(isreal) douR[iBELOW] = stkR[iBELOW];
		else       douC[iBELOW] = stkC[iBELOW];
	    }

	}

	if(isreal) sf_floatwrite  (douR, nBELOW, out);
	else       sf_complexwrite(douC, nBELOW, out);
    
    }

    if(isreal){
	free(dinR);
	free(douR);
	free(stkR);
    } else {
	free(dinC);
	free(douC);
	free(stkC);
    }

    exit (0);
}
