/* 
Threshold float/complex inputs given a constant/varying
threshold level.

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

void thrsample(sf_file in, sf_file out, bool complex_data, char* mode, 
	       float thr);

int main(int argc, char* argv[])
{
  int n1, n2;
  bool complex_data = false;/* true if input data is complex */
  char* mode;
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
  if (!sf_histint(in,"n1",&n1)) 
    sf_error("No n1= in input");
  /* number of traces */
  n2 = sf_leftsize(in,1);
  
  /* read threshold level from the command line */
  float thr;
  bool thrflag;
  bool fthrflag;
  sf_file fthr = NULL;
  int n1fthr, n2fthr;
  
  thrflag = sf_getfloat("thr",&thr);
  /* threshold level (>0) */
  fthrflag = sf_getstring("fthr");
  /* varying threshold level (>0) */
  
  if (!fthrflag){ /* constant threshold level */
    if (!thrflag)
      sf_error("Need threshold level");
  }
  else{ /* varying threshold level */
    fthr = sf_input("fthr");
    if (SF_FLOAT != sf_gettype(fthr))
      sf_error("Need float varying threshold level");
    if (!sf_histint(fthr,"n1",&n1fthr)) 
      sf_error("No n1= in fthr");
    n2fthr = sf_leftsize(fthr,1);
    if ( n1 != n1fthr || n2 != n2fthr)
      sf_error("Size mismatch: fthr & input must be of same size");
    if (!thrflag)
      thr = 1;
  }

  if (thr<=0) sf_error("Threshold must be >0");
  
  mode = sf_getstring("mode");
  /* 'soft', 'hard', 'nng' (default: soft)*/
  if (mode == NULL) mode = "soft";
  
  /* loop over samples */
  float *fthrsample;
  float thrvalue = thr;
  fthrsample = sf_floatalloc (1);
  
  for (int i=0; i<n1*n2; i++){
    if (fthrflag){
      sf_floatread(fthrsample,1,fthr);
      if (*fthrsample<0)
	sf_error("Varying threshold level must be >=0");
      else
	thrvalue = (*fthrsample)*thr;
    }
    thrsample(in,out,complex_data,mode,thrvalue);
  }
 
  sf_close(); 
  exit(0);
}

void thrsample(sf_file in, sf_file out, bool complex_data, char* mode, 
	       float thr)
{
  float *isample=NULL;
  float *osample=NULL;
  float complex *icsample=NULL;
  float complex *ocsample=NULL;
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
    absample = sqrt(pow(creal(*icsample),2) + pow(cimag(*icsample),2));
    
    if (absample > thr){
      if (strcmp(mode,"soft") == 0)
	*ocsample =(*icsample/absample)*(absample-thr);
      if (strcmp(mode,"hard") == 0)
	*ocsample =*icsample;
      if (strcmp(mode,"nng") == 0)
	*ocsample = *icsample-(pow(thr,2)/ *icsample);
    }
    else
      *ocsample = 0.;
    
    sf_complexwrite(ocsample,1,out);
  }
  else{
    sf_floatread(isample,1,in);

    if (fabs(*isample)>thr){
      if (strcmp(mode,"soft") == 0){
	if (*isample>thr)
	  *osample = *isample-thr;
	else
	  *osample = *isample+thr;
      }
      if (strcmp(mode,"hard") == 0)
	*osample = *isample;
      if (strcmp(mode,"nng") == 0)
	*osample = *isample-(pow(thr,2)/ *isample);
    }
    else
      *osample=0;
    
    sf_floatwrite(osample,1,out);
  }
}
