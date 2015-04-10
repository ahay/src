/* Read Trace And Header (tah) from standard input and FILTER 

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

The sftahfilter program is designed to apply a filter. Trace and 
header data (tah) are read from from standard input (usually a pipe).
A filter is read from the command line or from a file.  If the filter
is read from a file, it can be a single filter, or one filter for each
trace.  Future enhancements would be to use trace headers to define
a location in the filter file and read that filter or even interpolate
a filter.  That would support shot or receiver dependent filter for
surface consistant decon.  Another enhancement would be to add 
frequency domain bandpass filters.  

EXAMPLE:

sftahsort input=shots-receivers-23900_headfix.rsf           \\
   sort="xline:600,601 offset"                              \\
| sftahfilter filt=dephase.rsf                              \\
| sftahmakeskey pkey=xline skey=cdpt                        \\
| sftahwrite                                                \\
  verbose=1                                                 \\
  label2=cdpt  o2=1 n2=100 d2=1                             \\
  label3=xline o3=600 n3=1 d3=1                             \\
  output=dephasecmps.rsf                                    \\
>/dev/null

sfgrey <mutecmps.rsf | sfpen

In this example a deterministic dephase filter is applied to a prestack
datafile.

The shot organized prestack cmp data, shots-receivers-23900_headfix.rsf 
are read in xline offset order by sftahsort program.  The headers are 
in the file shots-receivers-23900_headfix_hdr.rsf, the headers 
parameter default.  The headers are merged with the trace amplitudes 
and the tah data sent down the pipe to apply a filter.

The sftahfilter program applies a filter contained in the dephase.rsf
file.
The program sftahmakeskey is used to create a secondary key used 
in the following sftahwrite to define the location to wrte the trace 
in the output file. Sftahmakeskey makes a secondary key (skey=cdpt) 
the count the traces starting in the a primary key gather (pkey=xline).
The input traces gathered by xline by sftahsort. Sftahmakeskey sets 
cdpt to 1 when the trace has a new xline.  If the trace has the same 
xline as the previous trace cdpt is incremented

Sftahwrite writes the the trace data to dephasecmp.rsf and the headers
are written to the file mutecmp_hdr.rsf.  The order of the data in the 
file is defined by the cdpt and xline trace headers, so the  data order
is (time,cmpt,xline).  Finally, the output volume is displayed using
sfgrey.

PARAMETERS

   floats filter=NULL

        A list of floats that is the filter to convolve on the input 
	traces.

   string filter_file=NULL

       Name of an rsf file that contains the filter(s)

   int filter_index_t0

        Index of time=0 on the filter

   int  verbose=1       
        
        flag to control amount of print
        0 terse, 1 informative, 2 chatty, 3 debug
	
*/

/*
   Program change history:
   date       Who             What
   12/12/2014 Karl Schleicher Original program.  Derived from Mtahstack.c
                              influences by sfconv, suconv, sufilter
*/

#include <string.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"
#include "convolve.h"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;
  int numfilter;
  float* filter=NULL;
  char* filter_file_name;
  sf_file filter_file=NULL;
  int filt_indx_t0;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  float* outtrace=NULL;
  int itrace=0;
  char **list_of_floats;
  float d1;
  float o1;


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
  if(NULL!=(list_of_floats=sf_getnstring("filter",&numfilter))){
    filter=sf_floatalloc(numfilter);
    if(!sf_getfloats("filter",filter,numfilter))
      sf_error("unable to read tmute");
  }

  filter_file_name = sf_getstring("filter_file");
  /* \n
     name of the rsf file containing the filter(s)
  */
  if(filter_file_name==NULL && filter==NULL){
    sf_error("your must input filter or filter_file");
  }
  if(filter_file_name!=NULL && filter!=NULL){
    sf_error("Define the filter once. Do not input filter AND filter_file.");
  }

  if(filter_file_name!=NULL){
    filter_file=sf_input("filter_file");
    if (!sf_histint(filter_file,"n1",&numfilter)){
      sf_error("No n1= in filter_file.");
    }
    filter=sf_floatalloc(numfilter);
    sf_floatread (filter,numfilter,filter_file);
  }

  if(!sf_getint("filt_indx_t0",&filt_indx_t0))filt_indx_t0=0;
    /* \n
     indx to time 0 in the filter.  Must be in the range [0,numfilter)
    */
  if(filt_indx_t0<0 || filt_indx_t0>=numfilter){
    fprintf(stderr,"length of the filter, numfilter=%d\n",numfilter);
    fprintf(stderr,"filt_indx_t0=%d numfilter=%d\n",filt_indx_t0,numfilter);
    sf_error("filt_indx_t0 muxt be in the range [0,numfilter)");
  }

  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);
  if(verbose>0){
    fprintf(stderr,"allocate outtrace.  length=%d\n",n1_traces+numfilter);
  }
  outtrace= sf_floatalloc(n1_traces+numfilter);


  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);

  /* this is how you get a header indx: indx_of_offset=segykey("offset"); */

  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data does not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data does not define o1");

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
 
  itrace=0;
  while (!(get_tah(intrace, fheader, n1_traces, n1_headers, in))){
    if(verbose>1 || (verbose==1 && itrace<5)){
      fprintf(stderr,"process tah %d in sftahmute\n",itrace);
    }

    /* this is how you read a filter for each location 
       if (!each) sf_floatread (ff,nf,filt);
    */
    convolve_init(numfilter,filter);
    convolve_lop (false, false,n1_traces,n1_traces+numfilter,intrace,outtrace);

    /********************/
    /* process the tah. */
    /********************/
    /* this program applies a mute the to the top of a trace */

    /* this is how you get a header value 
    if(typehead == SF_INT)offset=((int  *)fheader)[indx_of_offset];
    else                  offset=((float*)fheader)[indx_of_offset];
    */
    
    put_tah(&outtrace[filt_indx_t0], fheader, n1_traces, n1_headers, out);
    itrace++;
  }

  exit(0);
}

  
