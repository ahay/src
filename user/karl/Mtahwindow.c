/* Trace And Header WINDOW

tah is the abbreviation of Trace And Header.  Madagascar programs 
that begin with sftah are a designed to:
1- read trace and headers from separate rsf files and write them to 
   standard output (ie sftahread)
2- filter programs that read and write standard input/output and 
   process the tah data (eg sftahnmo, sftahstack)
3- read tah data from standard input and write separate rsf files for 
   the trace and headers data (ie sftahwrite)

These programs allow Seismic Unix (su) like processing in Madagascar.  
Some programs have su like names.

Some programs in this suite are sftahread, sftahgethw, ftahhdrmath, 
and sftahwrite.

The sftahwindow program select a subset of the traces.  Currently it will select
a subset of the time range of each trace.  In the future it needs to be 
upgraded to select traces based on trace headers.  It is modelled on suwind.

EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_field.rsf \\
| sftahwindow tmax=4.092 \\
   verbose=0  \\
| sftahwrite  \\
  verbose=1 \\
  output=npr3_field1.rsf \\
  mode=seq \\
>/dev/null

The headers are in the file npr3_field_hdr.rsf, the headers parameter 
default.  The headers are merged with the trace amplitudes and the tah 
data sent down the pipe for sftahwindow.  The trace is shortened to 
2047 samples to remove the two bad amplitudes observed at the end of 
most of the traces on this file.  The traces are sent to STDOUT to 
sftahwrite, which write the data sequentially to the output file (ie 
the output files is just a bunch of traces.

PARAMETERS
   Float tmax= maximum time in the input trace amplitude file.
        Maximum time in seconds to limit the output trace 

*/

/*
  Copyright (C) 2013 University of Texas at Austin
  
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
/*
   Program change history:
   date       Who             What
   09/15/2014 Karl Schleicher Original program
   06/10/2016 Karl Schleicher change time window parameter to tmax
*/
#include <string.h>
#include <limits.h>
#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  float d1;
  float o1;
  float tmax;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  int* iheader=NULL;
  float* intrace=NULL;
  int ns_out;
  int indx_ns; /*trace header index of ns key */
  char* keyname;
  int indx_key=0; /*trace header index of key identified by user's key parm */
  int key_min; /*minimum for header 'user key parameter to output */
  int key_max; /*maximum for header 'user key parameter to output' */
  int keyvalue; /* value of the header key on this trace */
  int *reject=NULL; /*list of bad user key headers*/
  int numreject;
  int ireject;
  bool rejecttrace;

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  sf_init (argc,argv);

  if(!sf_getint("verbose",&verbose))verbose=0;
      /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  sf_warning("verbose=%d",verbose);
 
  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/

  if(verbose>0)fprintf(stderr,"read in file name\n");  
  in = sf_input ("in");

  if(verbose>0)fprintf(stderr,"read out file name\n");
  out = sf_output ("out");

  if (!sf_histint(in,"n1_traces",&n1_traces))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n1_headers",&n1_headers)) 
    sf_error("input data does not define n1_headers");

  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data not define o1");

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  fheader = sf_floatalloc(n1_headers);
  iheader = (int*)fheader;
    
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  /* read and parse the user parameters */
  if(verbose>0)fprintf(stderr,"get the parameter tmax\n");
  if(!sf_getfloat("tmax",&tmax))tmax=(n1_traces-1)*d1+o1;
  /* maximum time in seconds to limit the output trace */
  ns_out=(tmax-o1)/d1+1.5;
  if(verbose>1){
    fprintf(stderr,"ns_out=%d, tmax=%f, o1=%f, d1=%f\n",ns_out, tmax, o1, d1);
  }

  if(ns_out>n1_traces){
    fprintf(stderr,
	    "input trace is shorter than computed from input paramater tmax\n");
    fprintf(stderr,"n1_traces=%d, ns_out=%d\n",n1_traces,ns_out);
    sf_error("input trace length is shorter than ns_out");
  }
  
  if(verbose>0)fprintf(stderr,"get the parameter reject\n");
  numreject=0;
  if(numreject>0){
    reject=sf_intalloc(numreject);
    sf_getints("reject",reject,numreject);
  }
  if (verbose>1){
    for (ireject=0; ireject<numreject; ireject++){
      fprintf(stderr,"reject[%d]=%d\n",ireject,reject[ireject]);
    }
  }

  if(verbose>0)fprintf(stderr,"get the parameters min, max, key\n");
  if(!sf_getint("min",&key_min))key_min=INT_MIN;
  if(!sf_getint("max",&key_max))key_max=INT_MAX;
  keyname=sf_getstring("key");
  fprintf(stderr,"key=%s\n",keyname);

  /**********************************************************/
  /* end code block for standard tah Trace And Header setup */
  /* continue with any sf_puthist this tah program calls to */
  /* add to the history file                                */
  /**********************************************************/

  sf_putint(out,"n1_traces",ns_out);

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);
  if(keyname!=NULL)indx_key=segykey(keyname);
  if(verbose>0)fprintf(stderr,"indx_key=%d\n",indx_key);
  /* kls upgrade this code to select based on trace header.
     loop for (i=0; i < n1_headers; i++) {
       look for segykeyword(i) in the input user parameters.  Count 
       numwindowheadernames
     }
     allocate array windowheadernames[numwindowheadernames], windowmin,windowmax
     loop and populate windowheadernames, windowmin, and windowmax.
     windowmin and windowmax should be double so large integers and
     floats in headers can be handled without rounding.  user input will
     look like tracf=1,1 or offset=0,600.

     as each trace is read, test header to determine if the traces should 
     be output or thrown away.

     update program to use tmax=3.500 instead of ns to apply time window.
     Karl Schleicher 05/22/2016
  */ 
  indx_ns=segykey("ns");


  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahgethw\n");
    /********************/
    /* process the tah. */
    /********************/
    rejecttrace=false;
    /* this program rejects traces if header key is not in range 
       [key_min,key_max] or is trase is in the reject list */
    if(keyname!=NULL){
      if(typehead == SF_INT){
	keyvalue=iheader[indx_key];
      } else {
	keyvalue=fheader[indx_key];
      }
      if(keyvalue<key_min || keyvalue>key_max)rejecttrace=true;
      
      for (ireject=0; ireject<numreject; ireject++){
	if(keyvalue==reject[ireject]){
	  rejecttrace=true;
	  break;
	}
      }
    
      if(verbose>1){
	  fprintf(stderr,"keyvalue=%d,key_min=%d,key_max=%d,",
		  keyvalue   ,key_min   ,key_max);
	  if(rejecttrace)fprintf(stderr,"rejecttrace=true\n");
	  else fprintf(stderr,"rejecttrace=false\n");
      }
    }

    if(!rejecttrace){
      if(verbose>2)fprintf(stderr,"ns_out=%d,indx_ns=%d\n",ns_out,indx_ns);
      if(typehead == SF_INT){
	iheader[indx_ns]=ns_out;
      } else {
	fheader[indx_ns]=(float)ns_out;
      }
      /***************************/
      /* write trace and headers */
      /***************************/
      put_tah(intrace, fheader, ns_out, n1_headers, out);
    }
  }

  exit(0);
}

  
