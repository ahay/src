/* Read data from an RSF file.
 *
 * MATLAB usage: rsf_write(data,file[,same])
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

#include <string.h>

#include "file.h"
#include "getpar.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int taglen, status, argc=2, i, ndim, len;
    const int *dim;
    size_t nbuf = BUFSIZ, nd, j;
    char *tag, *argv[] = {"matlab","-"}, *par, *filename;
    double *dr;
    float *p;
    char buf[BUFSIZ], key[5];
    bool same;
    FILE *file2;
    sf_file file;

    /* Check for proper number of arguments. */
    if (nrhs < 2 || nrhs > 3) mexErrMsgTxt("Two or three inputs required.");
    
    /* First input must be a string. */
    if (!mxIsChar(prhs[1]))
	mexErrMsgTxt("First input must be a string.");

    /* First input must be a row vector. */
    if (mxGetM(prhs[1]) != 1)
	mexErrMsgTxt("First input must be a row vector.");
    
    /* Get the length of the input string. */
    taglen = mxGetN(prhs[1]) + 1;

    /* Allocate memory for input string. */
    tag = mxCalloc(taglen, sizeof(char));

    /* Copy the string data from prhs[1] into a C string. */
    status = mxGetString(prhs[1], tag, taglen);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    if (3 == nrhs) {
        /* Input 3 must be a string. */
	if (!mxIsChar(prhs[2]))
	    mexErrMsgTxt("Input 3 must be a string.");

	/* Input 3 must be a row vector. */
	if (mxGetM(prhs[2]) != 1)
	    mexErrMsgTxt("Input 3 must be a row vector.");
	
	/* Get the length of the input string. */
	len = mxGetN(prhs[2]) + 1;
	
	/* Allocate memory for input string. */
	par = mxCalloc(len, sizeof(char));

	/* Copy the string data from prhs[2] into a C string. */
	status = mxGetString(prhs[2], par, len);
	if (status != 0) 
	    mexWarnMsgTxt("Not enough space. String is truncated.");

	same = (0 == (strncmp(par,"same",4)));
    } else {
	same = false;
    }

    sf_init(argc,argv);

    if (same) {
	file = sf_input(tag);
	filename = sf_histstring(file,"in");
	sf_fileclose(file);

	if (NULL == filename) mexErrMsgTxt("No in= in file.");
	file2 = fopen(filename,"ab");
	if (NULL == file2) 
	    mexErrMsgTxt("Could not open binary file for writing.");
    } else {
	file = sf_output(tag);
	file2 = NULL;
    }

    /* Input 2 must be a number. */
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Input 2 must be double.");
    
    /* data pointers */
    dr = mxGetPr(prhs[0]);
    
    /* get data dimensions */
    ndim=mxGetNumberOfDimensions(prhs[0]);
    dim=mxGetDimensions(prhs[0]);

    /* get data size */
    nd = mxGetNumberOfElements(prhs[0]);

    if (!same) {
	sf_setformat(file,mxIsComplex(prhs[0])?"native_complex":"native_float");
   
	/* Output */
	for (i=0; i < ndim; i++) {
	    sprintf(key,"n%d",i+1);
	    sf_putint(file,key,dim[i]);
	}
    }
	
    p = (float*) buf;

    for (j=0, nbuf /= sizeof(float); nd > 0; nd -= nbuf) {
	if (nbuf > nd) nbuf=nd;

	for (i=0; i < nbuf; i++, j++) {
	    p[i] = (float) dr[j];
	}

	if (same) {
	    if (nbuf != fwrite(p,sizeof(float),nbuf,file2)) 
		mexWarnMsgTxt("Writing problems.");
	} else {
	    sf_floatwrite(p,nbuf,file);
	}
    }

    if (same) {
	fclose(file2);
    } else {
	sf_fileclose(file);
    }
}
