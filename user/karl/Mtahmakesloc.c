/* Trace And Header MAKE SLOC KEY.

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

The sftahmakeslockey program will make a sloc keywhic is  useful for 
surface consistant scaling, decon, and residal statics.  An sloc is a 
surface location, ie a shot or a receiver location.  Programs could 
be written to drive off the group location (gx,gy), but it is better 
to just have an integer group surface location counts the group 
location from 1 to number group locations.  THat allows you to make 
shot record displays where the horizontal coordinate the group 
location index (a handy display to look for consistant receiver statics 
on shots).  You specify input surface coodinate, slocxy=gx,gy and output 
header key, slockey=tracf.

EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_gathers.rsf \\
| sftahmakeslockey slocxy=gx,gy slockey=tracf verbose=1 \\
| sftahwrite \\
   verbose=1                         \\
   mode=seq \\
   output=npr3_tracf.rsf \\
>/dev/null

In this example the cmp sorted prestack data, npr3_gathers.rsf,  are 
read by sftahread.  The headers are in the file npr3_gathers_hdr.rsf, 
the headers parameter default.  The headers are merged with the trace 
amplitudes and the tah data sent down the pipe for sftahmakeskey.
sftahmakeskey creates the cdpt header and sftahwrite creates a 4 
dimensional file.

PARAMETERS
   strings slocxy=

        list of 2 header keys that are the input surface location key
        so use to compute the slockey.

   string slockey

        the name of the output sloc header.
*/

/*
   Program change history:
   date       Who             What
   09/21/2015 Karl Schleicher Original program.  Derived from Mtahmakeskey.c
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
  char* slockey;
  int indx_tracex;
  int indx_tracey;
  int indx_sloc;
  double this_x;
  double this_y;
  double *tracex;
  double *tracey;
  int *slocarray;
  int slocarray_indx;

  int size_tracexy_array;
  int num_tracexy;

  int sloc;
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
 
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  if(verbose>0)fprintf(stderr,"call list of keys\n");
 
  list_of_keys=sf_getnstring("slocxy",&numkeys);
  /* two keys that are the surface x,y location to use to compute the trace
     slockey */

  if(list_of_keys==NULL)
    sf_error("The required parameter \"slocxy\" was not found.");
  /* I wanted to use sf_getstrings, but it seems to want a colon seperated
     list of keys (eg key=offset:ep:fldr:cdp) and I wanted a comma seperated
     list of keys (eg key=offset:ep:fldr:cdp).
  */
  if(numkeys!=2) sf_error("There must be exactly two header names in slocxy");

  /* print the list of keys */
  if(verbose>=1){
    fprintf(stderr,"numkeys=%d\n",numkeys);
    fprintf(stderr,"input trace xcoord=%s\n",list_of_keys[0]);
    fprintf(stderr,"input trace ycoord=%s\n",list_of_keys[1]);

    for(ikey=0; ikey<numkeys; ikey++){
      fprintf(stderr,"list_of_keys[%d]=%s\n",ikey,list_of_keys[ikey]);
    }
  }
  if(NULL==(slockey=sf_getstring("slockey")))
    sf_error("the required parameter \"slockey\" was not found");
  /* The name of the sloc key created by the program. */
  
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
  indx_tracex=segykey(list_of_keys[0]);
  indx_tracey=segykey(list_of_keys[1]);
  indx_sloc=segykey(slockey);

  size_tracexy_array=100;
  num_tracexy=0;

  tracex=malloc(size_tracexy_array*sizeof(double));
  tracey=malloc(size_tracexy_array*sizeof(double));
  slocarray=malloc(size_tracexy_array*sizeof(int));

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
 
  itrace=0;
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahgethw\n");
    /********************/
    /* process the tah. */
    /********************/
    /* this program computes a secondary key.  If one of the headers in 
       pkey changes, skey is set to 1.  Otherwise skey is the previous skey+1
    */

    if(typehead == SF_INT){
      this_x=((int*)fheader    )[indx_tracex];
      this_y=((int*)fheader    )[indx_tracey];
    } else {
      this_x=       fheader     [indx_tracex];
      this_y=       fheader     [indx_tracey];
    }
    slocarray_indx=tahinsert_unique_xy(&tracex,&tracey,&slocarray,
		                       &num_tracexy,&size_tracexy_array,
		                       this_x, this_y);
    sloc=slocarray[slocarray_indx];
    
    if(typehead == SF_INT) ((int*)fheader)[indx_sloc]=sloc;
    else                          fheader [indx_sloc]=sloc;
    
    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(intrace, fheader, n1_traces, n1_headers, out);
    itrace++;
  }

  exit(0);
}

  
