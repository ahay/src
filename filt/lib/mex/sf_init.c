#include <mex.h>

#include <rsf.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    char *argv[] = {"matlab","-"};
    int argc = 2;

    sf_init(argc,argv);
}
