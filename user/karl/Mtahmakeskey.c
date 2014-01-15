/* Trace And Header MAKE Secondary KEY.

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

The sftahmakeskey program will make a secondary key.  You sometimes
need a secondary key to number the traces in a gather.  For example 
after sorting the data in iline,xline,offset you usually cannot
store the data using the offset key because the offset sampling is 
irregular.  sftahmakeskey can be used to build the cdpt header from 
the iline and xline keys.  The secondary key, skey, will start with 
1 when a new iline,xline is encounterred.  As long as iline,xline
does not change, he skey will increase by 1 on each successive trace
until the iline,xline changes.  The output data can be stored in a 
file indexed by cdpt,xline,iline.

EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_gathers.rsf \\
| sftahmakeskey key=iline,xline skey=cdpt verbose=1 \\
| sftahwrite \\
   verbose=1                         \\
   label2="cdpt"  o2=1 n2=24  n2=1   \\
   label3="xline" o3=1 n3=188 d3=1   \\
   label4="iline" o4=1 n4=10 d4=1   \\
   output=mappedgather.rsf \\
>/dev/null

sfgrey <mappedgather.rsf | sfpen

In this example the cmp sorted prestack data, npr3_gathers.rsf,  are 
read by sftahread.  The headers are in the file npr3_gathers_hdr.rsf, 
the headers parameter default.  The headers are merged with the trace 
amplitudes and the tah data sent down the pipe for sftahmakeskey.
sftahmakeskey creates the cdpt header and sftahwrite creates a 4 
dimensional file.

PARAMETERS
   strings key= no default

        list of header keys to monitor to determine when to break 
	between gathers.  A gather is a sequence of traces with the 
	same value for all the header keys.  Stack summs traces in 
	the gather, divides by the fold, and outputs the stack trace.c
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
  float* fprevheader=NULL;
  int numkeys;
  int ikey;
  char** list_of_keys;
  int *indx_of_keys;
  char* skey;
  int indx_of_skey;
  int skeyvalue;
  bool pkeychanged;
  int itrace=0;

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  sf_init (argc,argv);

  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
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
  fprevheader = sf_floatalloc(n1_headers);
 
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  if(verbose>0)fprintf(stderr,"call list of keys\n");
 
  list_of_keys=sf_getnstring("pkey",&numkeys);
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
  if(NULL==(skey=sf_getstring("skey")))
    sf_error("the required parameter \"skey\" was not found");
  /* The name of the secondary key created by the program. */
  
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
  indx_of_skey=segykey(skey);

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
  skeyvalue=0;
 
  itrace=0;
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahgethw\n");
    /********************/
    /* process the tah. */
    /********************/
    /* this program computes a secondary key.  If one of the headers in 
       pkey changes, skey is set to 1.  Otherwise skey is the previous skey+1
    */
    pkeychanged=false;
    if(itrace>0){
      for(ikey=0; ikey<numkeys; ikey++){
	if(typehead == SF_INT){
	  if(((int*)fheader    )[indx_of_keys[ikey]]!=
	     ((int*)fprevheader)[indx_of_keys[ikey]]){
	    pkeychanged=true;
	    break;
	  }
	} else {
	  if(fheader[indx_of_keys[ikey]]!=fprevheader[indx_of_keys[ikey]]){
	    pkeychanged=true;
	    break;
	  }
	}
      }
    }
    if(pkeychanged) skeyvalue=1;
    else skeyvalue++;
    if(typehead == SF_INT) ((int*)fheader)[indx_of_skey]=skeyvalue;
    else                          fheader [indx_of_skey]=skeyvalue;
    
    if(skeyvalue==1){
      /* this is a new pkey, save the header so you know when it changes */
      memcpy(fprevheader,fheader,n1_headers*sizeof(int));
    }
    
    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(intrace, fheader, n1_traces, n1_headers, out);
    itrace++;
  }

  exit(0);
}

  
