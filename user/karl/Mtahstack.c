/* Read Trace And Header (tah) from STDIN, stack by mathing header key

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

The sftahstack program is designed to stack sorted data. Trace and 
header data (tah) are read from from standard input (usually a pipe).
Gathers to be stacked are defined by the input parameter \'key\'.
The key parameter is a list of header keys to watch.  A new gather 
starts when one of the header keys change.  On each gather the stack
trace is initialized to zero, traces in the gather are summed, and
the stack trace is divided by the time variant fold, and the stacked
trace written to output.  The output trace header is taken from the
first trace in the input gather.    

EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_gathers.rsf \\
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
amplitudes and the tah data sent down the pipe for nmo and stack.  The
sftahstack program uses both the iline and xline keys to determine
which traces blong to a gather.  Using both keys avoids a problem on 
edges of a survey when uising xline along may merge gathers across 
ilines (a special case that does sometimes happen). sftahwrite writes
the trace data to mappedstack.rsf and the headers are written to the
file mappedstack_hdr.rsf.  The order of the data in the output file
is defined by the iline and xline trace headers, so the  data order
is (time,xline,iline).  Finally, the output volume is displayed using
sfgrey.

PARAMETERS
   strings key= no default

        list of header keys to monitor to determine when to break 
	between gathers.  A gather is a sequence of traces with the 
	same value for all the header keys.  Stack summs traces in 
	the gather, divides by the fold, and outputs the stack trace.
*/

/*
   Program change history:
   date       Who             What
   11/15/2012 Karl Schleicher Original program.  Derived from Mtahgethw.c
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
  float* intrace=NULL;
  int numkeys;
  int ikey;
  char** list_of_keys;
  int *indx_of_keys;
  int indx_time;
  bool pkeychanged;
  float* stktrace=NULL;
  float* stkheader=NULL;
  float* time_variant_fold=NULL;
  bool eof_get_tah;
  int fold;
  int itrace=0;

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

  if(verbose>0)fprintf(stderr,"call list of keys\n");
 
  if(!(list_of_keys=sf_getnstring("key",&numkeys)))
    fprintf(stderr,"key not input\n");
  /*  List of header keys to monitor to determine when to break between 
      gathers.  A gather is a sequence of traces with the same value for
      all the header keys.  Stack summs traces in the gather, divides
      by the fold, and outputs the stack trace. */

  if(verbose>0)fprintf(stderr,"after sf_getnstring\n");

  if(list_of_keys==NULL)
    sf_error("The required parameter \"key\" was not found.");
  /* I wanted to use sf_getstrings, but it seems to want a colon seperated
     list of keys (eg key=offset:ep:fldr:cdp) and I wanted a comma seperated
     list of keys (eg key=offset:ep:fldr:cdp).
  numkeys=sf_getnumpars("key");
  if(numkeys==0)
    sf_error("The required parameter \"key\" was not found.");
  fprintf(stderr,"alloc list_of_keys numkeys=%d\n",numkeys);
  list_of_keys=(char**)sf_alloc(numkeys,sizeof(char*)); 
  sf_getstrings("key",list_of_keys,numkeys);
  */
  /* print the list of keys */
  if(verbose>1){
    fprintf(stderr,"numkeys=%d\n",numkeys);
    for(ikey=0; ikey<numkeys; ikey++){
      fprintf(stderr,"list_of_keys[%d]=%s\n",ikey,list_of_keys[ikey]);
    }
  }
  
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
  indx_of_keys=sf_intalloc(numkeys);
  for (ikey=0; ikey<numkeys; ikey++){
    /* kls need to check each of these key names are in the segy header and
       make error message for invalid keys.  Of does segykey do this? NO, just
       segmentation fault. */
    indx_of_keys[ikey]=segykey(list_of_keys[ikey]);
  }
  stktrace          = sf_floatalloc(n1_traces);
  stkheader         = sf_floatalloc(n1_headers);
  time_variant_fold = sf_floatalloc(n1_traces);

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
 
  itrace=0;
  eof_get_tah=get_tah(intrace, fheader, n1_traces, n1_headers, in);
  fold=0;
  while (!eof_get_tah){
    if(verbose>1 && itrace<5)fprintf(stderr,"process the tah in sftahstack\n");
    /********************/
    /* process the tah. */
    /********************/
    /* this program stacks sequential traces with matching pkey.  If one of the 
       headers in pkey changes, the summed trace is divided by the time 
       variant fold and output.
    */
    if(fold==0){
      memcpy(stkheader,fheader,n1_headers*sizeof(int));
      memcpy(stktrace ,intrace,n1_traces*sizeof(int));
      for(indx_time=0; indx_time<n1_traces; indx_time++){
	if(intrace[indx_time]!=0)time_variant_fold[indx_time]=1.0;
	else time_variant_fold[indx_time]=0.0;
      }
    } else {
      for(indx_time=0; indx_time<n1_traces; indx_time++){
	if(intrace[indx_time]!=0){
	  stktrace[indx_time]+=intrace[indx_time];
	  time_variant_fold[indx_time]++;
	}
      }
    }
    fold++;
    eof_get_tah=get_tah(intrace, fheader, n1_traces, n1_headers, in);

    /* did any of the header keys in indx_of_keys change? */
    pkeychanged=false;
    if(itrace>0){
      for(ikey=0; ikey<numkeys; ikey++){
	if(typehead == SF_INT){
	  if(((int*)fheader  )[indx_of_keys[ikey]]!=
	     ((int*)stkheader)[indx_of_keys[ikey]]){
	    pkeychanged=true;
	    break;
	  }
	} else {
	  if(fheader[indx_of_keys[ikey]]!=stkheader[indx_of_keys[ikey]]){
	    pkeychanged=true;
	    break;
	  }
	}
      }
    }
    /* if one of the headers changes, apply fold recovery, output trace, and 
       set fold=0.  Fold=0 will initialize the stktrace to zero at top f loop*/
    if(pkeychanged){
      /***********************************/
      /* divide by the time variant fold */
      /***********************************/
      if(verbose>1)fprintf(stderr,"pkeychanged.  divide by fold\n");
      for(indx_time=0; indx_time<n1_traces; indx_time++){
	if(time_variant_fold[indx_time]<1.0)time_variant_fold[indx_time]=1.0;
	stktrace[indx_time]/=time_variant_fold[indx_time];
      }   
      /* kls set any headers? Maybe fold? offset? */
      /***************************/
      /* write trace and headers */
      /***************************/
      if(verbose>1)fprintf(stderr,"put_tah\n");
      put_tah(stktrace, stkheader, n1_traces, n1_headers, out);
      fold=0;
    }
    itrace++;
  }

  exit(0);
}

  
