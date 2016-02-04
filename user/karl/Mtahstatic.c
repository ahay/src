/*  Trace And Header STATIC

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

The sftahstatic applies a static computed using the sstat and gstat trace
headers. These headers are the source and receiver static in milliseconds.
The program also uses the input parameter sign. The static applied:

tin=tout+sign*(sstat+gstat)

This means positive sstat+gstat and sign=1 will shift a trace "up".

Interpolation error is less than 1% for frequencies less than 60% of
the Nyquist frequency. 
                                                                             
EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_gathers.rsf \\
| sftahstatic sign=-1 \\
| sftahnmo \\
   verbose=1  \\
   tnmo=0,.373,.619,.826,.909,1.017,1.132,1.222,1.716,3.010 \\
   vnmo=9086,10244,11085,10803,10969,11578,12252,12669,14590,17116 \\
| sftahstack key=iline,xline verbose=1 \\
| sftahwrite \\
   verbose=1                           \\
   label2="xline" o2=1 n2=188 d2=1   \\
   label3="iline" o3=1 n3=345 d3=1   \\
   output=mappedstack.rsf \\
>/dev/null

sfgrey <mappedstack.rsf | sfpen

In this example the cmp sorted prestack data, npr3_gathers.rsf,  are 
read by sftahread.  The headers are in the file npr3_gathers_hdr.rsf, 
the headers parameter default.  The headers are merged with the trace 
amplitudes and the tah data sent down the pipe for statics, nmo, and 
stack.  The sftahstack program uses both the iline and xline keys to 
determine which traces blong to a gather.  Using both keys avoids a 
problem on edges of a survey when uising xline along may merge gathers 
across ilines (a special case that does sometimes happen). sftahwrite 
writes the trace data to mappedstack.rsf and the headers are written 
to the file mappedstack_hdr.rsf.  The order of the data in the output 
file is defined by the iline and xline trace headers, so the  data 
order is (time,xline,iline).  Finally, the output volume is displayed 
using sfgrey.
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
   08/26/2015 Karl Schleicher Original program based on Mtahnmo.c
   02/04/2016 Karl Schleicher Write selfdoc
*/
#include <string.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>
#include <math.h>

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
  float* intrace=NULL;
  float o1;
  float d1;
  int indx_sstat;
  int indx_gstat;
  int indx_time;
  float* tin_of_tout=NULL;
  float* outtrace=NULL;
  int itrace;
  float tstat;
  int static_sign;

  /******************************************************/
  /* code block for standard tah Trace And Header setup */
  /******************************************************/

  sf_init (argc,argv);

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
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
 
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  /* maybe I should add some validation that n1== n1_traces+n1_headers+2
     and the record length read in the second word is consistent with 
     n1.  */
  /**********************************************************/
  /* end code block for standard tah Trace And Header setup */
  /* continue with any sf_puthist this tah program calls to */
  /* add to the history file                                */
  /**********************************************************/

  if(!sf_getint("sign",&static_sign))static_sign=1;

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/
  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data does not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data does not define o1");
  /* Kls should read label1 and verify it is time  
  if (!sf_histstring(in,"label1",&label1))
    sf_error("input data not define label1");
  */

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);
  /* get index to keys I will be using */
  indx_sstat=segykey("sstat");
  indx_gstat=segykey("gstat");
  /* kls what other header keys do I use?  inline? xline? cdp? */

  if(verbose>0)fprintf(stderr,"allocate arrays for the trace loop\n");
  tin_of_tout=sf_floatalloc(n1_traces);
  outtrace         =sf_floatalloc(n1_traces);
  

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahnmo\n");
    /********************/
    /* process the tah. */
    /********************/
    /* this program applies moveout */
    /* kls read the statics from the headers */
    if(typehead == SF_INT){   
      /* just cast the header to int  */
      tstat=((int*)fheader)[indx_sstat]+((int*)fheader)[indx_gstat];
    } else {
      tstat=       fheader [indx_sstat]+       fheader [indx_gstat];
    }
    
    /* kls su code is:
       for (itime=0; itime<nt; ++itime)
       t[itime] = tmin + itime*dt + sign*tstat;
       // args are:
       // nt input, dt input, tmin in, indata
       // value to use for points before input tmin, after tmax, nt out, 
       //            output timed, output trace 
       ints8r(nt, dt, tmin, intrace.data,
              0.0, 0.0, nt, t, outtrace.data);
    */	

    /* it would be computationally faster to convolve the interp filter 
       to the trace, but this is easy copy of sustatic.  Cost of ints8r
       should no be too bad and I want to get initial code up quickly */
    itrace=0;
    for(indx_time=0; indx_time<n1_traces; indx_time++){
      tin_of_tout[indx_time]=indx_time*d1+static_sign*tstat/1000;
      if(itrace==0 && verbose>4){
	fprintf(stderr,"tin_of_tout[%d]=%f\n",
		indx_time   , tin_of_tout[indx_time]);
      }
    }
    ints8r(n1_traces,d1,0,
	   intrace,0.0,0.0,
	   n1_traces,tin_of_tout,
	   outtrace);

    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(outtrace, fheader, n1_traces, n1_headers, out);
    itrace++;
  }
  
  exit(0);
}

  
