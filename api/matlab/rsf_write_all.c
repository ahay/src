/* Write complete RSF file, both header and data, in one call.
 *
 * MATLAB usage:
 *   rsf_write_all(file,cmdargs,data[,delta[,origin[,label[,unit]]]])
 *
 *   Note! Only cmdargs={'out=stdout'} was tested
 *
 * Written by Henryk Modzelewski, UBC EOS SLIM
 * Created February 2012
 */
/*
  Copyright (C) 2012 The University of British Columbia at Vancouver
  
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

#include <stdio.h>
#include <string.h>
#include <mex.h>
#include <rsf.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int strlen, status, argc, i, odim, Odim;
    mwSize ndim;
    const mwSize *dim=NULL;
    size_t nbuf = BUFSIZ, nd, j;
    char *strtag=NULL, **argv, *filename=NULL;
    double *dr=NULL, *di=NULL;
    double *ddlt=NULL, *dorg=NULL;
    mxArray *pca;
    float *p=NULL;
    sf_complex *c=NULL;
    char buf[BUFSIZ], key[8];
    sf_file file=NULL;
    
    /* Check for proper number of arguments. */
    if (nlhs==0 && nrhs==0) {
	printf("RSF_WRITE_ALL Writes complete RSF file from MATLAB\n");
	printf("Usage:\n");
	printf("\trsf_write_all(file,cmdargs,data[,delta[,origin[,label[,unit]]]])\n");
	printf("Where:\n");
	printf("\tInput:\n");
	printf("\t\tfile is the RSF-file name\n");
	printf("\t\tcmdargs holds cell array of strings with RSF comman-line options\n");
	printf("\t\t\t(only cmdargs={'out=stdout'} was tested)\n");
	printf("\t\tdata holds the data array\n");
	printf("\tOptional output:\n");
	printf("\t\tdelta holds row vector of d# float values for header\n");
	printf("\t\t\tlength(delta) must be >= ndims(data)\n");
	printf("\t\torigin holds row vector of o# float values for header\n");
	printf("\t\tlabel holds cell array of label# string for from header\n");
	printf("\t\tunit holds cell array of unit# string values for header\n");
	printf("\t\t\tlength(origin|label|unit) must be = length(delta)\n");
	printf("\tNote: the file is closed on exit from this function\n");
    	return;
    }
    if (nrhs < 3 || nrhs > 7)
	 mexErrMsgTxt("3 to 7 inputs required:\n\tfile,cmdargs,data[,delta[,origin[,label[,unit]]]]\n\t Note! Only cmdargs={'out=stdout'} was tested.");
    if (nlhs != 0)
	 mexErrMsgTxt("This function has no outputs");

    /* File name must be a string. */
    if (!(mxIsChar(prhs[0])&&mxGetM(prhs[0])==1)) mexErrMsgTxt("File name must be a string.");
    /* Get the length of the input string. */
    strlen = mxGetN(prhs[0]) + 1;
    /* Allocate memory for input string. */
    filename = mxCalloc(strlen, sizeof(char));
    /* Copy the string data from prhs[0] into a C string. */
    status = mxGetString(prhs[0], filename, strlen);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    /* Command line options must be a cell array of strings */
    if (!mxIsCell(prhs[1])) mexErrMsgTxt("Labels must be a cell array.");
    odim = mxGetNumberOfElements(prhs[1]);
    argc = 2;
    argv = mxCalloc((argc+odim),sizeof(char*));
    argv[0] = mxCalloc(7,sizeof(char)); sprintf(argv[0],"matlab");
    argv[1] = mxCalloc(2,sizeof(char)); sprintf(argv[1],"-");
    for (i=0; i<odim; i++) {
    	pca = mxGetCell(prhs[1], i);
	if (!(mxIsChar(pca)&&mxGetM(pca)==1)) mexErrMsgTxt("Command-line option must be a string.");
	strlen = mxGetN(pca) + 1;
	argv[argc] = mxCalloc(strlen, sizeof(char));
	status = mxGetString(pca, argv[argc], strlen);
	argc++;
    }

    /* Data must be a double array. */
    if (!mxIsDouble(prhs[2])) mexErrMsgTxt("Data must be a double array.");
    /* get data dimensions */
    ndim=mxGetNumberOfDimensions(prhs[2]);
    dim=mxGetDimensions(prhs[2]);
    /* set default # of output dimansions for metadata */
    Odim = ndim;
    /* get data size */
    nd = mxGetNumberOfElements(prhs[2]);

    /* Deltas must be a double row vector. */
    if (nrhs > 3) {
	if (!mxIsDouble(prhs[3])) mexErrMsgTxt("Delta must be double.");
	if (mxGetM(prhs[3]) != 1) mexErrMsgTxt("Deltas must be a row vector.");
	Odim = mxGetN(prhs[3]);
	ddlt = mxGetPr(prhs[3]);
	if (Odim < ndim) mexErrMsgTxt("Deltas does not have enough elements.");
    }

    /* Origins must be a double row vector. */
    if (nrhs > 4) {
	if (!mxIsDouble(prhs[4])) mexErrMsgTxt("Origin must be double.");
	if (mxGetM(prhs[4]) != 1) mexErrMsgTxt("Origins must be a row vector.");
	odim = mxGetN(prhs[4]);
	dorg = mxGetPr(prhs[4]);
	if (odim != Odim) mexErrMsgTxt("Origins has wrong number of elements.");
    }

    /* Labels must be a cell array of strings. */
    if (nrhs > 5) {
	if (!mxIsCell(prhs[5])) mexErrMsgTxt("Labels must be a cell array.");
	odim = mxGetNumberOfElements(prhs[5]);
	if (odim != Odim) mexErrMsgTxt("Labels has wrong number of elements.");
    }

    /* Units must be a cell array of strings. */
    if (nrhs > 6) {
	if (!mxIsCell(prhs[6])) mexErrMsgTxt("Units must be a cell array.");
	odim = mxGetNumberOfElements(prhs[6]);
	if (odim != Odim) mexErrMsgTxt("Units has wrong number of elements.");
    }

    sf_init(argc,argv);
    file = sf_output(filename);
    sf_setformat(file,mxIsComplex(prhs[2])?"native_complex":"native_float");

    /* Write header for at least ndim dimensions */
    for (i=0; i < Odim; i++) {
	/* sizes */
        sprintf(key,"n%d",i+1);
        if (i<ndim) sf_putint(file,key,(int)dim[i]);
        else sf_putint(file,key,1);
	/* deltas */
        sprintf(key,"d%d",i+1);
        if (nrhs > 3) sf_putfloat(file,key,(float)ddlt[i]);
        else sf_putfloat(file,key,1.);
	/* origins */
        sprintf(key,"o%d",i+1);
        if (nrhs > 4) sf_putfloat(file,key,(float)dorg[i]);
        else sf_putfloat(file,key,0.);
	/* labels */
	if (nrhs > 5) {
	    pca = mxGetCell(prhs[5], i);
	    if (!mxIsChar(pca)) mexErrMsgTxt("Label must be a string.");
	    strlen = mxGetN(pca) + 1;
	    strtag = mxCalloc(strlen, sizeof(char));
	    status = mxGetString(pca, strtag, strlen);
	    sprintf(key,"label%d",i+1);
	    if (strlen > 1) sf_putstring(file,key,strtag);
	    mxFree(strtag);
	}
	/* units */
	if (nrhs > 6) {
	    pca = mxGetCell(prhs[6], i);
	    if (!mxIsChar(pca)) mexErrMsgTxt("Unit must be a string.");
	    strlen = mxGetN(pca) + 1;
	    strtag = mxCalloc(strlen, sizeof(char));
	    status = mxGetString(pca, strtag, strlen);
	    sprintf(key,"unit%d",i+1);
	    if (strlen > 1) sf_putstring(file,key,strtag);
	    mxFree(strtag);
	}
    }
    
    /* Write data */
    if (mxIsComplex(prhs[2])) {
	/* complex data */
	c = (sf_complex*) buf;

	dr = mxGetPr(prhs[2]);

	/* pointer to imaginary part */
	di = mxGetPi(prhs[2]);
	
	for (j=0, nbuf /= sizeof(sf_complex); nd > 0; nd -= nbuf) {
	    if (nbuf > nd) nbuf=nd;
	    
	    for (i=0; i < nbuf; i++, j++) {
		c[i] = sf_cmplx((float) dr[j],(float) di[j]);
	    }
	    
	    sf_complexwrite(c,nbuf,file);
	}

    } else { 
	/* real data */
	p = (float*) buf;

	dr = mxGetPr(prhs[2]);

	for (j=0, nbuf /= sizeof(float); nd > 0; nd -= nbuf) {
	    if (nbuf > nd) nbuf=nd;
	    
	    for (i=0; i < nbuf; i++, j++) {
		p[i] = (float) dr[j];
	    }
	    
	    sf_floatwrite(p,nbuf,file);
	}
    } 

    sf_fileclose(file);
    sf_close();

    for (i=0;i<argc;i++) mxFree(argv[i]); mxFree(argv);
}
