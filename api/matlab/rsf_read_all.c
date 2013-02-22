/* Read complete RSF file, both header and data, in one call.
 *
 * MATLAB usage:
 *   [data [size [delta [origin [label [unit]]]]]] = rsf_read_all(file)
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
#include <sys/stat.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
    int strlen, status, argc=2, dim, n[SF_MAX_DIM], i, esize; /* len; */
    size_t nbuf = BUFSIZ, nd, j;
    char *strtag, *argv[] = {"matlab","-"}; /* *par; */
    double *pr, *pi=NULL;
    char buf[BUFSIZ];
/*    static off_t shift=0; */
    sf_datatype type;
    sf_file file;
    char key[8];
    int nn[SF_MAX_DIM];
    struct stat fstat;

    /* Check for proper number of arguments. */
    if (nlhs==0 && nrhs==0) {
	printf("RSF_READ_ALL Reads complete RSF file into MATLAB\n");
	printf("Usage:\n");
	printf("\t[data [size [delta [origin [label [unit]]]]]] = rsf_read_all(file)\n");
	printf("Where:\n");
	printf("\tInput:\n");
	printf("\t\tfile is the RSF-file name\n");
	printf("\tOutput:\n");
	printf("\t\tdata holds the data array\n");
	printf("\tOptional output:\n");
	printf("\t\tsize is size(data); vector of n# int values from header\n");
	printf("\t\tdelta holds vector of d# float values from header\n");
	printf("\t\torigin holds vector of o# float values from header\n");
	printf("\t\tlabel holds cell array of label# string values from header\n");
	printf("\t\tunit holds cell array of unit# string values from header\n");
	printf("\tNote: the file is closed on exit from this function\n");
    	return;
    }
    if (nrhs != 1) mexErrMsgTxt("1 input required: file");
    if (nlhs < 0 || nlhs > 6)
	 mexErrMsgTxt("1 to 6 outputs required:\n\tdata[,size[,delta[,origin[,label[,unit]]]]]");

    /* File name must be a string. */
    if (!(mxIsChar(prhs[0])&&mxGetM(prhs[0])==1)) mexErrMsgTxt("Second input must be a string.");
    /* Get the length of the input string. */
    strlen = mxGetN(prhs[0]) + 1;
    /* Allocate memory for input string. */
    strtag = mxCalloc(strlen, sizeof(char));
    /* Copy the input filename into a C string. */
    status = mxGetString(prhs[0], strtag, strlen);
    if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

    sf_init(argc,argv);
    if (stat(strtag, &fstat))
       {printf("FATAL ERROR: file %s does not exist\n",strtag); sf_close(); return;}
    file = sf_input(strtag);
    type = sf_gettype (file);
    esize = sf_esize(file);
    dim = sf_filedims(file,n);
    for (i=0; i < dim; i++) {
	sprintf(key,"n%d",i+1);
     	if (!sf_histint(file,key,&nn[i]))
	   {printf("FATAL ERROR: missing %s in file %s\n",key,strtag); sf_close(); return;}
    }
    for (i=dim; i < SF_MAX_DIM; i++) {
	sprintf(key,"n%d",i+1);
	if (sf_histint(file,key,&nn[i]) && nn[i]>0) dim++;
	else break;
    }

    /* Crete data array and pointers */
    plhs[0] = mxCreateNumericArray( dim, nn, mxDOUBLE_CLASS,
        type==SF_COMPLEX?mxCOMPLEX:mxREAL);
    pr = mxGetPr(plhs[0]);
    if (type==SF_COMPLEX) pi = mxGetPi(plhs[0]);
    nd = mxGetNumberOfElements(plhs[0]);

    /* Return size */
    if (nlhs > 1) {
	double *p;
    	plhs[1] = mxCreateDoubleMatrix(1,dim,mxREAL);;
	p = mxGetPr(plhs[1]);
	for (i=0; i < dim; i++) p[i] = nn[i];
    }

    /* Return deltas */
    if (nlhs > 2) {
	float dummy; double *p;
    	plhs[2] = mxCreateDoubleMatrix(1,dim,mxREAL);;
	p = mxGetPr(plhs[2]);
	for (i=0; i < dim; i++) {
	    sprintf(key,"d%d",i+1);
     	    sf_histfloat(file,key,&dummy);
	    p[i] = dummy;
	}
    }

    /* Return origins */
    if (nlhs > 3) {
	float dummy; double *p;
    	plhs[3] = mxCreateDoubleMatrix(1,dim,mxREAL);;
	p = mxGetPr(plhs[3]);
	for (i=0; i < dim; i++) {
	    sprintf(key,"o%d",i+1);
     	    sf_histfloat(file,key,&dummy);
	    p[i] = dummy;
	}
    }

    /* Return labels */
    if (nlhs > 3) {
    	plhs[4] = mxCreateCellMatrix(1,dim);
	for (i=0; i < dim; i++) {
	    sprintf(key,"label%d",i+1);
	    mxSetCell(plhs[4], i,
	        mxCreateString(sf_histstring(file,key)));
	}
    }

    /* Return units */
    if (nlhs > 4) {
    	plhs[5] = mxCreateCellMatrix(1,dim);
	for (i=0; i < dim; i++) {
	    sprintf(key,"unit%d",i+1);
	    mxSetCell(plhs[5], i,
	        mxCreateString(sf_histstring(file,key)));
	}
    }

    /* Read data */
    for (j=0, nbuf /= esize; nd > 0; nd -= nbuf) {
	if (nbuf > nd) nbuf=nd;

	switch(type) {
	    case SF_FLOAT:

		sf_floatread((float*) buf,nbuf,file);
		for (i=0; i < nbuf; i++, j++) {
		    pr[j] = (double) ((float*) buf)[i];
		}
		break;
	    case SF_INT:

		sf_intread((int*) buf,nbuf,file);
		for (i=0; i < nbuf; i++, j++) {
		    pr[j] = (double) ((int*) buf)[i];
		}
		break;
	    case SF_COMPLEX:

		sf_complexread((sf_complex*) buf,nbuf,file);
		for (i=0; i < nbuf; i++, j++) {
		    pr[j] = (double) crealf(((sf_complex*) buf)[i]);
		    pi[j] = (double) cimagf(((sf_complex*) buf)[i]);
		}
		break;
	    default:
		mexErrMsgTxt("Unsupported file type.");
		break;
	}
    }

    sf_fileclose(file);
    sf_close();
}
