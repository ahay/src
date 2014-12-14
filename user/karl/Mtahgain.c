/* Read Trace And Header (tah) from standard input and apply GAIN

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

The sftahgain program is designed to apply various amplitude gain 
functions.  These gain functions include multiply by a constant, 
trace balance, clipping, AGC, t**pow, exp(epow*t) and AGC.  Trace 
and header data (tah) are read from standard input (usually a pipe).
Data is scaled:
out(t) = scale * BAL{CLIP[AGC{[t^tpow * exp(epow * t)}]}
Then trace and header data is written to standard output.

EXAMPLE:

sftahsort input=shots-receivers-23900_headfix.rsf          \\
   sort="xline:600,601 offset"                             \\
| sftahgain agc=.5                                         \\
| sftahmakeskey pkey=xline skey=cdpt                       \\
| sftahwrite                                               \\
  verbose=1                                                \\
  label2=cdpt  o2=1 n2=100 d2=1                            \\
  label3=xline o3=600 n3=1 d3=1                            \\
  output=agccmps.rsf                                       \\
>/dev/null

sfgrey <agccmps.rsf | sfpen

In this example the shot organized prestack data in the file 
shots-receivers-23900_headfix.rsf are read in xline offset order by 
sftahsort program.  The headers are in the file 
shots-receivers-23900_headfix_hdr.rsf, the headers parameter default.
The headers are merged with the trace amplitudes and the tah data sent 
down the pipe for scaling (agc).

sftahagc applies a .5 second automatic gain control (agc)

The program sftahmakeskey is used to create a secondary key used 
in the following sftahwrite to define the location to wrte the trace 
in the output file. Sftahmakeskey makes a secondary key (skey=cdpt) 
the count the traces starting in the a primary key gather (pkey=xline).
The input traces gathered by xline by sftahsort. Sftahmakeskey sets 
cdpt to 1 when the trace has a new xline.  If the trace has the same 
xline as the previous trace cdpt is incremented

Sftahwrite writes the the trace data to agccmp.rsf and the headers are 
written to the file agccmp_hdr.rsf.  The order of the data in the output 
file is defined by the cdpt and xline trace headers, so the  data order
is (time,cmpt,xline).  Finally, the output volume is displayed using
sfgrey.

PARAMETERS
   int verbose=1
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug

   float scale=1.0
        multiply data by this float

   float tpow=0.0
        multiply data by t^tpow

   float epow=0.0
        multiply data by exp(epow*t)

   float agc=0.0
        Length of agc window in seconds.  0.0 means no agc

   Operation order:

   out(t) = scale * BAL{CLIP[AGC{[t^tpow * exp(epow * t)}]}	
*/

/*
   Program change history:
   date       Who             What
   12/12/2014 Karl Schleicher Original program.  Derived from Mtahmute.c
                              with lots of inspiration from sugain.
*/
#include <string.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"

int main(int argc, char* argv[])
{
  int verbose;

  float scale;
  float tpow;
  float epow;
  float agc;
  int agc_samples;

  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  int indx_time;
  int itrace=0;
  int numxstart;
  int numtstart;
  char **list_of_floats;
  float* xstart;
  float* tstart;
  int indx_of_offset;
  float offset;
  float d1;
  float o1;
  float* tpow_epow;
  float agc_start;
  int iagc_start;


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
 
  if(!sf_getfloat("scale",&scale))scale=1.0;
  /* multiply data by this float */

  if(!sf_getfloat("tpow",&tpow))tpow=0.0;
  /* multiply data by t^tpow */

  if(!sf_getfloat("epow",&epow))epow=0.0;
  /* multiply data by exp(epow*t) */

  if(!sf_getfloat("agc",&agc))agc=0.0;
  /* Length of agc window in seconds.  0.0 means no agc */

  numxstart=0;
  numtstart=0;

  if(NULL!=(list_of_floats=sf_getnstring("xstart",&numxstart))){
    /* \n
       list of offsets that correspont to tstart and define the offset 
       dependent start time for the agc scaling */
    xstart=sf_floatalloc(numxstart);
    if(!sf_getfloats("xxstart",xstart,numxstart))
      sf_error("unable to read xstart");
  }
  if(NULL!=(list_of_floats=sf_getnstring("tstart",&numtstart))){
    /* \n
       list of times that correspont to xstart and define the offset 
       dependent start time for the agc scaling */
    tstart=sf_floatalloc(numtstart);
    if(!sf_getfloats("tstart",tstart,numtstart))sf_error("unable to read tstart");
  }
  if(numxstart!=numtstart)sf_error("bad start parameters: numxstart!=numtstart");

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

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);

  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data does not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data does not define o1");

  /* initialize scaling arrays that do not change */

  /* operation is:
     out(t) = scale * BAL{CLIP[AGC{[t^tpow * exp(epow * t) * 
                                   ( in(t)-bias )]^gpow}
                              ]
                         }
  */
  tpow_epow=sf_floatalloc(n1_headers);
  if(tpow==0){
    for (indx_time=0; indx_time<n1_traces; indx_time++){
      tpow_epow[indx_time]=1.0;
    }
  } else {
    for (indx_time=0; indx_time<n1_traces; indx_time++){
      tpow_epow[indx_time]=pow(o1+indx_time*d1,tpow);
    }
  }
  if(epow!=0){
    for (indx_time=0; indx_time<n1_traces; indx_time++){
      tpow_epow[indx_time]*=exp(epow*(o1+indx_time*d1));
    }
  }

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
 
  itrace=0;
  while (!(get_tah(intrace, fheader, n1_traces, n1_headers, in))){
    if(verbose>1 || (verbose==1 && itrace<5)){
      fprintf(stderr,"process tah %d in sftahmute\n",itrace);
    }
    /********************/
    /* process the tah. */
    /********************/
    /* this program applies a start the to the top of a trace */

    /* apply tpow and epow */
    if(tpow!=0 || epow!=0){
      for(indx_time=0; indx_time<n1_traces;
	  indx_time++){
	intrace[indx_time]*=tpow_epow[indx_time];
      }
    }
    /* apply agc */
    if(agc!=0){
      sf_error("agc option for sftahgain still needs to be writen");
      /* this program applies a start the to the top of a trace */
      
      if(typehead == SF_INT)offset=((int  *)fheader)[indx_of_offset];
      else                  offset=((float*)fheader)[indx_of_offset];
      intlin(numxstart,xstart,tstart,
	     tstart[0],tstart[numxstart-1],1,
	     &offset,&agc_start);
      if(agc_start<o1)agc_start=o1;
      
      iagc_start=(int)(((agc_start-o1)/d1)+.5);
      if(0)fprintf(stderr,"istart_start=%d\n",iagc_start);
      
    }

    put_tah(intrace, fheader, n1_traces, n1_headers, out);
    itrace++;
  }

  exit(0);
}

  
