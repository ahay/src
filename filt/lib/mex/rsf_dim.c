#include <mex.h>

#include <rsf.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int taglen, status, argc=2, dim, n[SF_MAX_DIM], i;
    char *tag, *argv[] = {"matlab","-"};
    double *p;
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

    sf_init(argc,argv);
    file = sf_input(tag);
    dim = sf_filedims(file,n);
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
    p = mxGetPr(plhs[0]);

    for (i=0; i < dim; i++) {
	p[i] = n[i];
    }
}
