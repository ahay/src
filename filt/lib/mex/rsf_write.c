/* Read data from an RSF file.
 *
 * MATLAB usage: rsf_write(file,data)
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
    int taglen, status, argc=2, i, ndim;
    const int *dim;
    size_t nbuf = BUFSIZ, nd, j;
    char *tag, *argv[] = {"matlab","-"};
    double *dr, *di;
    float *p;
    char buf[BUFSIZ], key[5];
    sf_file file;

    /* Check for proper number of arguments. */
    if (nrhs != 2) mexErrMsgTxt("Two inputs required.");
    
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
    file = sf_output(tag);

    /* Input 2 must be a number. */
    if (!mxIsDouble(prhs[1])) mexErrMsgTxt("Input 2 must be double.");

    /* data pointers */
    dr = mxGetPr(prhs[1]);
    di = mxGetPr(prhs[1]);
    
    /* get data dimensions */
    ndim=mxGetNumberOfDimensions(prhs[1]);
    dim=mxGetDimensions(prhs[1]);

    /* get data size */
    nd = mxGetNumberOfElements(prhs[1]);

    sf_setformat(file,mxIsComplex(prhs[1])?"native_complex":"native_float");
   
    /* Output */
    for (i=0; i < ndim; i++) {
	sprintf(key,"n%d",i+1);
	sf_putint(file,key,dim[i]);
    }

    p = (float*) buf;

    for (j=0, nbuf /= sizeof(float); nd > 0; nd -= nbuf) {
	if (nbuf > nd) nbuf=nd;

	for (i=0; i < nbuf; i++, j++) {
	    p[i] = (float) dr[j];
	}

	sf_floatwrite(p,nbuf,file);
    }

    sf_fileclose(file);
}
