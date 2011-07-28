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

static void thrsample(sf_file in, sf_file out, bool complex_data, 
		      const char* mode, float thr);

int main(int argc, char* argv[])
{
    off_t n1, n2, i, n1fthr, n2fthr;
    float thr, thrvalue;
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
    fthrflag = (bool) (NULL != sf_getstring("fthr"));
    /* varying threshold level (positive number) */
    
    if (!fthrflag){ /* constant threshold level */
	if (!thrflag)
	    sf_error("Need threshold level");
    }
    else{ /* varying threshold level */
	fthr = sf_input("fthr");
	if (SF_FLOAT != sf_gettype(fthr))
	    sf_error("Need float varying threshold level");
	if (!sf_histlargeint(fthr,"n1",&n1fthr)) 
	    sf_error("No n1= in fthr");
	n2fthr = sf_leftsize(fthr,1);
	if ( n1 != n1fthr || n2 != n2fthr)
	    sf_error("Size mismatch: fthr & input must be of same size");
	if (!thrflag)
	    thr = 1;
    }
    
    if (thr<0) sf_error("Threshold must be >=0");
    
    mode = sf_getstring("mode");
    /* 'soft', 'hard', 'nng' (default: soft)*/
    if (mode == NULL) mode = "soft";
    
    /* loop over samples */
    thrvalue = thr;
    fthrsample = sf_floatalloc (1);
    
    for (i=0; i<n1*n2; i++){
	if (fthrflag){
	    sf_floatread(fthrsample,1,fthr);
	    if (*fthrsample<0)
		sf_error("Varying threshold level must be >=0");
	    else
		thrvalue = (*fthrsample)*thr;
	}
	thrsample(in,out,complex_data,mode,thrvalue);
    }
    
 
    exit(0);
}

static void thrsample(sf_file in, sf_file out, bool complex_data, 
		      const char* mode, float thr)
{
    float *isample=NULL;
    float *osample=NULL;
    sf_complex *icsample=NULL;
    sf_complex *ocsample=NULL;
    float absample;
    
    if (complex_data){
	icsample = sf_complexalloc(1);
	ocsample = sf_complexalloc(1);
    }
    else{
	isample = sf_floatalloc(1);
	osample = sf_floatalloc(1);
    }
    
    if (complex_data){
	sf_complexread(icsample,1,in);
	absample = cabsf(*icsample);
	
	if (absample > thr){
	    if (strcmp(mode,"soft") == 0) {
#ifdef SF_HAS_COMPLEX_H
		*ocsample =(*icsample/absample)*(absample-thr);
#else
		*ocsample = sf_crmul(*icsample,(absample-thr)/absample);
#endif
		
	    } else if (strcmp(mode,"hard") == 0) {
		*ocsample =*icsample;
	    } else if (strcmp(mode,"nng") == 0) {
#ifdef SF_HAS_COMPLEX_H
		*ocsample = *icsample-(thr*thr/ *icsample);
#else
		*ocsample = sf_csub(*icsample,
				    sf_cdiv(sf_cmplx(thr*thr,0.),*icsample));
#endif
	    }
	} else {
	    *ocsample = sf_cmplx(0.,0.);
	}
	
	sf_complexwrite(ocsample,1,out);
    }
    else{
	sf_floatread(isample,1,in);
	
	if (fabsf(*isample)>thr){
	    if (strcmp(mode,"soft") == 0){
		if (*isample>thr)
		    *osample = *isample-thr;
		else
		    *osample = *isample+thr;
	    }
	    if (strcmp(mode,"hard") == 0)
		*osample = *isample;
	    if (strcmp(mode,"nng") == 0)
		*osample = *isample-(thr*thr/ *isample);
	}
	else
	    *osample=0;
	
	sf_floatwrite(osample,1,out);
    }
    free(isample);
    free(osample);
    free(icsample);
    free(ocsample);
}
