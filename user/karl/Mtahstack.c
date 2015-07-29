/* Read Trace And Header (tah) from STDIN, stack matching header keys

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
    tnmo=0,2,6,10.5,16 vnmo=1500,1500,2250,3250,3700  \\
| sftahstack key=iline,xline \
    xmute=0,20000 tmute=0,20 ntaper=25 \\
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

   floats xmute= NULL

        List of floats the same length as list of floats in the tmute
	parameter.  The (xmute,tmute) pairs are interpolated using the
	trace headers offset to determine trace start time.  The mute is
	always increased to the first non-zero sample.  The default mutes 
	at the first non-zero sample.

   floats tmute= NULL

        List of floats the same length as list of floats in the xmute
	parameter.  The (xmute,tmute) pairs are interpolated using the
	trace headers offset to determine trace start time. The mute is
	always increased to the first non-zero sample.  The default mutes 
	at the first non-zero sample.

   float ntaper=12
        the length of the taper to use at the start of the trace.
	
*/

/*
   Program change history:
   date       Who             What
   11/15/2012 Karl Schleicher Original program.  Derived from Mtahgethw.c
*/
#include <string.h>
#include <rsf.h>
# include <rsf_su.h>
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
  int eof_get_tah;
  int fold;
  int itrace=0;
  int ntaper;
  int numxmute;
  int numtmute;
  float* taper;
  char **list_of_floats;
  float* xmute;
  float* tmute;
  int indx_of_offset;
  float offset;
  float d1;
  float o1;
  float mute_start;
  int imute_start;
  int indx_taper;


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
  /* get the mute parameters */
  if(NULL==(list_of_floats=sf_getnstring("xmute",&numxmute))){
    xmute=NULL;
    numxmute=0;
  } else {
    xmute=sf_floatalloc(numxmute);
    if(!sf_getfloats("xmute",xmute,numxmute))sf_error("unable to read xmute");
  }
  if(NULL==(list_of_floats=sf_getnstring("tmute",&numtmute))){
    tmute=NULL;
    numtmute=0;
  } else {
    tmute=sf_floatalloc(numtmute);
    if(!sf_getfloats("tmute",tmute,numtmute))sf_error("unable to read tmute");
  }
  if(numxmute!=numtmute)sf_error("bad mute parameters: numxmute!=numtmute");
  if(numxmute<=0) {
      ntaper=0;
      indx_of_offset=0;
      taper=NULL;
  } else {
    if(!sf_getint("ntaper",&ntaper))ntaper=12;
    /* \n
       length of the taper on the stack mute
    */
    taper=sf_floatalloc(ntaper);
    for(indx_time=0; indx_time<ntaper; indx_time++){
      float val_sin=sin((indx_time+1)*SF_PI/(2*ntaper));
      taper[indx_time]=val_sin*val_sin;
    }
    indx_of_offset=segykey("offset");
  }

  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data does not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data does not define o1");

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
      for(indx_time=0; indx_time<n1_traces; indx_time++){
	time_variant_fold[indx_time]=0.0;
	stktrace[indx_time]=0.0;
      }
    }
    if(xmute==NULL)mute_start=o1;
    else{
      if(typehead == SF_INT)offset=((int  *)fheader)[indx_of_offset];
      else                  offset=((float*)fheader)[indx_of_offset];
      intlin(numxmute,xmute,tmute,
	     tmute[0],tmute[numxmute-1],1,
	     &offset,&mute_start);
      if(mute_start<o1)mute_start=o1;
    }
    imute_start=(int)(((mute_start-o1)/d1)+.5);
    if(0)fprintf(stderr,"imute_start=%d\n",imute_start);
    for(; imute_start<n1_traces && intrace[imute_start]==0;
	  imute_start++); /* increate imute_start to first non-zero sample */
    if(0)fprintf(stderr,"updated imute_start=%d\n",imute_start);
    for(indx_time=imute_start, indx_taper=0; 
	indx_time<imute_start+ntaper && indx_time<n1_traces; 
	indx_time++, indx_taper++){
      stktrace[indx_time]+=taper[indx_taper]*intrace[indx_time];
      time_variant_fold[indx_time]+=taper[indx_taper];
    }
    for(; indx_time<n1_traces; indx_time++){
      stktrace[indx_time]+=intrace[indx_time];
      time_variant_fold[indx_time]++;
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

  
