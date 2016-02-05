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
| sftahwindow ns=2047 \\
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
   int ns= n1 in the input trace amplitude file.

        The number of samples in the input traces.

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
*/
#include <string.h>
#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  int* iheader=NULL;
  float* intrace=NULL;
  int ns;
  int indx_ns;
  
  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  sf_init (argc,argv);

  if(!sf_getint("verbose",&verbose))verbose=1;
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

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  fheader = sf_floatalloc(n1_headers);
  iheader = (int*)fheader;
    
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  if(verbose>0)fprintf(stderr,"get the parameter nt\n");
  if(!sf_getint("ns",&ns))ns=n1_traces;
  
  if(ns>n1_traces){
    fprintf(stderr,"input trace length is shorter than input paramter nt\n");
    fprintf(stderr,"n1_traces=%d, ns=%d\n",n1_traces,ns);
    sf_error("input trace length is shorter than input paramter nt");
  }
  
  /**********************************************************/
  /* end code block for standard tah Trace And Header setup */
  /* continue with any sf_puthist this tah program calls to */
  /* add to the history file                                */
  /**********************************************************/

  sf_putint(out,"n1_traces",ns);

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);
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
    /* this program prints selected header keys */
    
    if(verbose>2)fprintf(stderr,"ns=%d,indx_ns=%d\n",ns,indx_ns);
    if(typehead == SF_INT){
      iheader[indx_ns]=ns;
    } else {
      fheader[indx_ns]=(float)ns;
    }
    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(intrace, fheader, ns, n1_headers, out);
  }

  exit(0);
}

  
