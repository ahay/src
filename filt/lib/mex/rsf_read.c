/* Read data from an RSF file.
 *
 * MATLAB usage: array = rsf_read(file,dims)
 *
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

#include <mex.h>

#include <rsf.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int taglen, status, argc=2, dim, n[SF_MAX_DIM], i, esize;
    size_t nbuf = BUFSIZ, nd, j;
    char *tag, *argv[] = {"matlab","-"};
    double *p;
    char buf[BUFSIZ];
    sf_datatype type;
    sf_file file;

    /* Check for proper number of arguments. */
    if (nrhs != 1 && nrhs != 2) {
	mexErrMsgTxt("One or two inputs required.");
    } else if (nlhs > 1) { 
	mexErrMsgTxt("Too many output arguments.");
    }

    /* First input must be a string. */
    if (!mxIsChar(prhs[0]))
	mexErrMsgTxt("First input must be a string.");

    /* First input must be a row vector. */
    if (mxGetM(prhs[0]) != 1)
	mexErrMsgTxt("First input must be a row vector.");
    
    /* Get the length of the input string. */
    taglen = mxGetN(prhs[0]) + 1;

    /* Allocate memory for input string. */
    tag = mxCalloc(taglen, sizeof(char));

    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[0], tag, taglen);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    sf_init(argc,argv);
    file = sf_input(tag);

    dim = sf_filedims(file,n);
    type = sf_gettype (file);

    if (!sf_histint(file,"esize",&esize)) esize=sizeof(float);

    /* Output */
    if (1 == nrhs) { /* read everything */
	nd = 1;
	for (i=0; i < dim; i++) {
	    nd *= n[i];
	}

	plhs[0] = mxCreateNumericArray(dim,n,mxDOUBLE_CLASS,mxREAL);
	p = mxGetPr(plhs[0]);

	for (j=0, nbuf /= esize; nd > 0; nd -= nbuf) {
	    if (nbuf > nd) nbuf=nd;

	    switch(type) {
		case SF_FLOAT:
		    sf_floatread((float*) buf,nbuf,file);
		    for (i=0; i < nbuf; i++, j++) {
			p[j] = (double) ((float*) buf)[i];
		    }
		    break;
		case SF_INT:
		    sf_intread((int*) buf,nbuf,file);
		    for (i=0; i < nbuf; i++, j++) {
			p[j] = (double) ((int*) buf)[i];
		    }
		    break;
		default:
		    mexErrMsgTxt("Unsupported file type.");
		    break;  
	    }
	}
    }
}
