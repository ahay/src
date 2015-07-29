/* Get a parameter from a RSF file.
 *
 * MATLAB usage: par = rsf_par(file,name,type,default)
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
    int len, status, argc=2, i;
    char *tag, *name, *type, *argv[] = {"matlab","-"};
    double *par, *def;
    float f;
    bool b;
    sf_file file=NULL;

    /* Check for proper number of arguments. */
    if (nrhs != 4) {
	mexErrMsgTxt("Four inputs required.");
    } else if (nlhs > 1) { 
	mexErrMsgTxt("Too many output arguments.");
    }

    /* Process inputs */

    /* Input 1 must be a string. */
    if (!mxIsChar(prhs[0]))
	mexErrMsgTxt("Input 1 must be a string.");

    /* Input 1 must be a row vector. */
    if (mxGetM(prhs[0]) != 1)
	mexErrMsgTxt("Input 1 must be a row vector.");

    /* Get the length of the input string. */
    len = mxGetN(prhs[0]) + 1;

    /* Allocate memory for input string. */
    tag = mxCalloc(len, sizeof(char));

    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[0], tag, len);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    /* Input 2 must be a string. */
    if (!mxIsChar(prhs[1]))
	mexErrMsgTxt("Input 2 must be a string.");

    /* Input 2 must be a row vector. */
    if (mxGetM(prhs[1]) != 1)
	mexErrMsgTxt("Input 2 must be a row vector.");

    /* Get the length of the input string. */
    len = mxGetN(prhs[1]) + 1;

    /* Allocate memory for input string. */
    name = mxCalloc(len, sizeof(char));

    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[1], name, len);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    /* Input 3 must be a string. */
    if (!mxIsChar(prhs[2]))
	mexErrMsgTxt("Input 3 must be a string.");

    /* Input 3 must be a row vector. */
    if (mxGetM(prhs[2]) != 1)
	mexErrMsgTxt("Input 3 must be a row vector.");

    /* Get the length of the input string. */
    len = mxGetN(prhs[2]) + 1;

    /* Allocate memory for input string. */
    type = mxCalloc(len, sizeof(char));

    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[2], type, len);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");
 
    /* Input 4 must be a number. */
    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]))
	mexErrMsgTxt("Input 4 must be a number.");

    /* Input 4 must be a scalar. */
    if (mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 1)
	mexErrMsgTxt("Input 4 must be a scalar.");

    def = mxGetPr(prhs[0]);

    /* Output */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    par = mxGetPr(plhs[0]);
    *par = *def;

    /* Processing */

    sf_init(argc,argv);
    file = sf_input(tag);

    switch (type[0]) {
	case 'i': case 'd': /* integer */
	    if (sf_histint(file,name,&i)) *par = i;
	    break;
	case 'l': case 'b': /* logical */
	    if (sf_histbool(file,name,&b)) *par = b? 1.: 0.;
	    break;
	case 'f': case 'r': /* floating-point */
	    if (sf_histfloat(file,name,&f)) *par = f;
	    break;
	default:
	    mexErrMsgTxt("Unknown type.");
	    break;
    }

}
