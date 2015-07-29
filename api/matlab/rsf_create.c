/* Write RSF header with desired info to disk
 *
 * MATLAB usage: rsf_create(new_file_name, old_file_name | dims)
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
    int taglen, taglen2, status, argc=2, i, ndim;
    char *tag, *tag2, *argv[] = {"matlab","-"};
    double *dr;
    char key[5];
    sf_file newfile=NULL, oldfile=NULL;

    /* Check for proper number of arguments. */
    if (nrhs != 2) mexErrMsgTxt("Two inputs are required.");

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
    newfile = sf_output(tag);

    if (mxIsChar(prhs[1])) {
	/* Second input must be a row vector. */
	if (mxGetM(prhs[1]) != 1)
	    mexErrMsgTxt("Second input must be a row vector.");

	/* Get the length of the input string. */
	taglen2 = mxGetN(prhs[1]) + 1;

	/* Allocate memory for input string. */
	tag2 = mxCalloc(taglen2, sizeof(char));

	/* Copy the string data from prhs[1] into a C string. */
	status = mxGetString(prhs[1], tag2, taglen2);
	if (status != 0) 
	    mexWarnMsgTxt("Not enough space. String is truncated.");
	
	oldfile = sf_input(tag2);
    } else {
	if (mxGetN(prhs[1]) != 1)
	    mexErrMsgTxt("Second input must be a column vector.");

	/* Second input must be a number. */
	if (!mxIsDouble(prhs[1])) mexErrMsgTxt("Second input must be double.");

	ndim = mxGetM(prhs[1]);
	dr = mxGetPr(prhs[1]);

	for (i=0; i < ndim; i++) {
	    sprintf(key,"n%d",i+1);
	    sf_putint(newfile,key,(int) dr[i]);
	}
    }

    sf_setformat(newfile,"native_float");
    sf_fileflush(newfile,oldfile); /* The actual writing to disk */

}
