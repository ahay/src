#include <mex.h>

#include <rsf.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int taglen, status;
    char *tag;
    sf_file file;

    /* Check for proper number of arguments. */
    if (nrhs != 1) {
	mexErrMsgTxt("One input required.");
    } else if (nlhs > 1) { 
	mexErrMsgTxt("Too many output arguments.");
    }

    /* Input must be a string. */
    if (!mxIsChar(prhs[0]))
	mexErrMsgTxt("Input must be a string.");

    /* Input must be a row vector. */
    if (mxGetM(prhs[0]) != 1)
	mexErrMsgTxt("Input must be a row vector.");
    
    /* Get the length of the input string. */
    taglen = mxGetN(prhs[0]) + 1;

    /* Allocate memory for input string. */
    tag = mxCalloc(taglen, sizeof(char));

    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[0], tag, taglen);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    /* Output */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    /* file = (sf_file*) mxGetPr(plhs[0]); */

    file = sf_input(tag);
}
