/* Threshold float/complex inputs given a constant/varying threshold level.

Methods available:
- soft
- hard
- non-negative Garrote (nng)

Written by: Gilles Hennenfent & Colin Russell, UBC
Created: February 2006
*/

/*
  Copyright (C) 2006 The University of British Columbia at Vancouver
  
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

#include <rsf.h>

#include "thr.c"

int main(int argc, char* argv[])
{
    off_t n1, n2, i, n1fthr, n2fthr;
    float thr, thrvalue, *input, *output;
    bool thrflag;
    bool fthrflag;
    float *fthrsample;
    sf_file fthr = NULL;
    bool complex_data = false;/* true if input data is complex */
    const char* mode;
    sf_file in, out; /* Input and output files */
  
    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");
    
    /* check that input is either float or complex */
    if (SF_COMPLEX == sf_gettype(in)) 
	complex_data = true;
    else if (SF_FLOAT != sf_gettype(in)) 
	sf_error("Need float or complex input");
    
    /* number of time samples/trace */
    if (!sf_histlargeint(in,"n1",&n1)) 
	sf_error("No n1= in input");
    /* number of traces */
    n2 = sf_leftsize(in,1);
    
    /* read threshold level from the command line */
    thrflag = sf_getfloat("thr",&thr);
    /* threshold level (positive number) */
    //fthrflag = (bool) (NULL != sf_getstring("fthr"));
    /* varying threshold level (positive number) */
    
    //if (!fthrflag){ /* constant threshold level */
	if (!thrflag) sf_error("Need threshold level");
    //}
    //else{ /* varying threshold level */
	//fthr = sf_input("fthr");
	//if (SF_FLOAT != sf_gettype(fthr))
	 //   sf_error("Need float varying threshold level");
	//if (!sf_histlargeint(fthr,"n1",&n1fthr)) 
	  //  sf_error("No n1= in fthr");
	//n2fthr = sf_leftsize(fthr,1);
	//if ( n1 != n1fthr || n2 != n2fthr)
	  //  sf_error("Size mismatch: fthr & input must be of same size");
	//if (!thrflag)
	  //  thr = 1;
    //}
    
    if (thr<0) sf_error("Threshold must be >=0");
    
    mode = sf_getstring("mode");
    /* 'soft', 'hard', 'nng' (default: soft)*/
    if (mode == NULL) mode = "soft";
    
    thrvalue = thr;
    
    input = sf_floatalloc (n2*n1);

    output = sf_floatalloc (n2*n1);

    sf_floatread(input,n2*n1,in);
    
    thr_init(complex_data /* if input is complex */,
             n2 /* number of traces */,
             n1 /* number of time samples */,
	     thr /* thresholding value */,
	     mode /*thresholding mode */);

    thr_lop(input,output);

    sf_floatwrite(output,n2*n1,out);
  
    exit(0);
}
