/* Read Trace And Header from seperate files, combine, write to pipe

tah is the abbreviation of Trace And Header.  Madagascar programs 
that begin with sftah are a group of programs designed to:
1- read trace and headers from separate rsf files and write them to 
   standard output (ie sftahread)
2- filter programs that read and write standard input/output and process 
   the tah data (eg sftahnmo, sftahstack)
3- read tah data from standard input and write separate rsf files for the
   trace and headers data (ie sftahwrite)

These programs allow Seismic Unix (su) like processing in Madagascar.  
Some programs have su like names.

Some programs in this suite are sf_tahread, sf_tahgethw, f_tahhdrmath, 
and sf_tahwrite.
*/
/*
   Program change history:
   date       Who             What
   04/26/2013 Karl Schleicher Original program
*/

#include <string.h>

#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

/* sparingly make some global variables. */
int verbose;

int main(int argc, char* argv[])
{
  sf_file infile=NULL, out=NULL, inheaders=NULL;
  int n1_traces;
  int n1_headers;

  int n_traces, n_headers; 
  sf_datatype typehead;
  sf_datatype typein;
  float* header=NULL;
  float* intrace=NULL;
  int i_trace;
  char* infile_filename=NULL;
  char* headers_filename=NULL;

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
  if(verbose>0)fprintf(stderr,"read name of input file name\n");
  
  infile_filename=sf_getstring("input");
  if(infile_filename==NULL) infile = sf_input ("in");
  else infile = sf_input (infile_filename);

  if(verbose>0)
    fprintf(stderr,"set up output file for tah - should be stdout\n");
  out = sf_output ("out");

  if(verbose>0)fprintf(stderr,"read name of input headers file\n");
  if(!(headers_filename=sf_getstring("headers"))) {
    sf_error("headers parameter not input");
  } else {
    inheaders = sf_input(headers_filename);
  }

  if (!sf_histint(infile,"n1",&n1_traces))
    sf_error("input file does not define n1");
  if (!sf_histint(inheaders,"n1",&n1_headers)) 
    sf_error("input headers file does not define n1");

  n_traces=sf_leftsize(infile,1);
  n_headers=sf_leftsize(inheaders,1);
  if(verbose>0)
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
  /*fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);*/
  header = sf_floatalloc(n1_headers);
 
  /*fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);*/

  intrace= sf_floatalloc(n1_traces);

  /*fprintf(stderr,"need to add 2 words.  Record type and record length\n"); */
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
  {
    int tempint;
    tempint=segykey("ep");
    fprintf(stderr,"tempint=%d\n",tempint);
  }
  /* put the history from the input file to the output */
  sf_fileflush(out,infile);

  if(verbose>0)fprintf(stderr,"start trace loop n_traces=%d\n",n_traces);
  for (i_trace=0; i_trace<n_traces; i_trace++){
    if(verbose>1 ||(verbose>0 && i_trace<5)){
      fprintf(stderr,"i_trace=%d\n",i_trace);
    }
    /**************************/
    /* read trace and headers */
    /**************************/
    sf_floatread(header,n1_headers,inheaders);
    sf_floatread(intrace,n1_traces,infile);

    /***************************/
    /* write trace and headers */
    /***************************/
    /* trace and header will be preceeded with words:
       1- type record: 4 charactors 'tah '.  This will support other
          type records like 'htah', hidden trace and header.
       2- the length of the length of the trace and header. */
    put_tah(intrace, header, n1_traces, n1_headers, out);
  }

  exit(0);
}

  
