/* 
   tahread: Trace And Header READ.

   tah is the abbreviation of Trace And Header.  It identifies a group of
   programs designed to:
   1- read trace and headers from separate rsf files and write them to 
      standard output
   2- filter programs that read and write standard input/output and process 
      the tah data
   3- read tah data from standard input and write separate rsf files for the
      trace and headers data

   These programs allow Seismic Unix (su) like processing in Madagascar.  
   Some programs have su like names.

   Some programs in this suite are sf_tahread, sf_tahgethw, f_tahhdrmath, 
   and sf_tahwrite.
*/
/*
   Program change history:
   date       Who             What
   04/26/2012 Karl Schleicher Original program
*/

#include <string.h>
#include <rsf.h>

/* I do not know how to include this header or link to the right library
   obviously this is a terrible cludge */
#include "/home/karl/RSFSRC/system/seismic/segy.h"
#include "/home/karl/RSFSRC/system/seismic/segy.c"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file infile=NULL, out=NULL, inheaders=NULL;
  int n1_traces;
  int n1_headers;

  int n_traces, n_headers; 
  sf_datatype typehead;
  sf_datatype typein;
  float* fheader=NULL;
  float* intrace=NULL;
  int* iheader=NULL;
  int i_trace;
  int tempint;

  sf_init (argc,argv);

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
  fprintf(stderr,"verbose=%d\n",verbose);
 
  /*****************************************/
  /* initialize the input and output files */
  /*****************************************/
  fprintf(stderr,"read in file name\n");
  
  if(NULL!=sf_getstring("input")) infile = sf_input ("input");
  else                            infile = sf_input ("in");

  fprintf(stderr,"read out file name\n");
  out = sf_output ("out");

  fprintf(stderr,"read name of input headers file\n");
  if(!(inheaders = sf_input("headers")))sf_error("headers parameter not input");
  
  if (!sf_histint(infile,"n1",&n1_traces))
    sf_error("input file does not define n1");
  if (!sf_histint(inheaders,"n1",&n1_headers)) 
    sf_error("input headers file does not define n1");

  n_traces=sf_leftsize(infile,1);
  n_headers=sf_leftsize(inheaders,1);
  fprintf(stderr,"n_traces=%d, n_headers=%d\n",n_traces,n_headers);

  if(n_traces!=n_headers){
    fprintf(stderr,"n_traces=%d, n_headers=%d\n",n_traces,n_headers);
    sf_error("number of traces and headers must be the same");
  }

  typehead = sf_gettype(inheaders);
  if (SF_INT != typehead && SF_FLOAT != typehead ) 
    sf_error("Need float or int input headers.");

  typein = sf_gettype(infile);
  if (SF_FLOAT != typein ) 
    sf_error("input file must contain floats.");
  
  /* must be float or int */
  fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  if (SF_INT == typehead) iheader = sf_intalloc(n1_headers);
  else                    fheader = sf_floatalloc(n1_headers);
 
  fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);

  intrace= sf_floatalloc(n1_traces);

  fprintf(stderr,"need to add 2 words.  Record type and record length\n");
  sf_putint(out,"n1",n1_traces+n1_headers+2);
  sf_putint(out,"n1_traces",n1_traces);

  /*************************************************************************/
  /* put info into history so later tah programs can figure out the headers*/
  /*************************************************************************/
  sf_putint(out,"n1_headers",n1_headers);
  if(SF_INT == typehead) sf_putstring(out,"header_format","native_int");
    else                 sf_putstring(out,"header_format","native_float");
  sf_putstring(out,"headers",sf_getstring("headers"));
  /* the list of any extra header keys */
  segy_init(n1_headers,inheaders);
  segy2hist(out,n1_headers);

  /* put the history from the input file to the output */
  sf_fileflush(out,infile);

  fprintf(stderr,"start trace loop\n");
  for (i_trace=0; i_trace<n_traces; i_trace++){
    if(i_trace<5)fprintf(stderr,"i_trace=%d\n",i_trace);

    /**************************/
    /* read trace and headers */
    /**************************/
    if(SF_INT == typehead) sf_intread  (iheader,n1_headers,inheaders);
    else                   sf_floatread(fheader,n1_headers,inheaders);
    sf_floatread(intrace,n1_traces,infile);

    /***************************/
    /* write trace and headers */
    /***************************/
    /* trace and header will be preceeded with words:
       1- type record: 4 charactors 'tah '.  This will support other
          type records like 'htah', hidden trace and header.
       2- the length of the length of the trace and header. */
    sf_charwrite("tah ",4,out);
    tempint=n1_traces+n1_headers;
    sf_intwrite(&tempint, 1, out);
    if(SF_INT == typehead) sf_intwrite  (iheader,n1_headers,out);
    else                   sf_floatwrite(fheader,n1_headers,out);
    sf_floatwrite(intrace,n1_traces,out);
  }

  exit(0);
}

  
