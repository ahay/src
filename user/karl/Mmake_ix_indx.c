/* MAKE Iline Xline INDX files for quick 3D data subsets (superbins)

These indexes are intended to be used by future sftahsort and sftah5dinterp.

EXAMPLE:
Runnning the programs:

sfmake_ix_indx           \\
   verbose=1             \\
   input=npr3_fielda.rsf \\
   indxdir=npr3_field    \\
   ilineinc=10           \\
   iline=iline           \\
   xline=xline           \\
   append=no             \\
>/dev/null 
sfmake_ix_indx           \\
   verbose=1             \\
   input=npr3_fieldb.rsf \\
   indxdir=npr3_field    \\
   ilineinc=10           \\
   iline=iline           \\
   xline=xline           \\
   append=yes            \\
>/dev/null 

Will create a set of files:
npr3_field/il0
npr3_field/il5
...
npr3_field/il350
and
npr3_field/filename_indx

The file npr3_field/il15 will contain the trace headers of the traces with 
trace header "iline" nearest to 15.  The file also has the file number and
the trace number, so you can locate the trace in the input files (either 
npr3_fielda.rsf or npr3_fieldb.rsf.  The trace headers in the file are all 
sorted by xline.  With this information you can quickly find all the traces
that are in a range of ilines and xlines.  This supports sftah5dinterp which
processes all the traces in a midpoint superbin that might be 800 meters by 800
meters (about 16 ilines and 32 xlines).  This is not a simple sort problem 
because sftah5dinterp processes data in overlapping bin (i.e. the 800 meter 
superbin centers moveup by 400 meters).  Overlapping superbins are supported by
allowing traces to be reread.

The program also allows sftahsort to read from multiple files. This is useful
on larger 3D projects where the input data is on multiple segy files.  
Previously I merged the files into one big file after running sfsegyread.  This
required an additional copy of all the data to be saved on disk.    
*/

/*
   Program change history:
   date       Who             What
   11/07/2014 Karl Schleicher Original program based in Mtahsort
*/
#include <string.h>

/* includes for realpath */
#include <stdio.h>
#include <sys/param.h>
#include <stdlib.h>

#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"
#include "extraalloc.h"

/* sparingly make some global variables. */
int verbose;
int numsortheadernames;

typedef struct {
  int filenumber;
  int tracenumber;
  float keyvalue[1];
} indxkey;

int qsort_compare(const void *a, const void *b){
  int ikey;
  indxkey* a_indxkey=*(indxkey**)a;
  indxkey* b_indxkey=*(indxkey**)b;

  for (ikey=0; ikey<numsortheadernames; ikey++){
    /* as soon as keys differ you can return the compare */
    if(a_indxkey->keyvalue[ikey] > b_indxkey->keyvalue[ikey])return 1;
    if(a_indxkey->keyvalue[ikey] < b_indxkey->keyvalue[ikey])return -1;
  }
  /* keys are all equal.  */
  return 0;
}

int rounddouble(double num){
  return num<0 ? (int)(num-.5) : (int)(num+.5);
}

int main(int argc, char* argv[])
{
  sf_file infile=NULL, out=NULL, inheaders=NULL;
  int n1_traces;
  int n1_headers;

  int n_traces, n_headers; 
  sf_datatype typehead;
  sf_datatype typein;
  float* header=NULL;
  char* infile_filename=NULL;
  char* headers_filename=NULL;
  /*
  indxkey** pntrs_to_indxkeys;
  char* sort;
  char** sortheadernames;
  char** sorttokens1;
  double* signsortname;
  double* sortmin;
  double* sortmax;
  double* sortinc;
  int n_traces_sort;
  int n_pntrs_to_indxkeys;
  int i1_headers;
  int* indx_of_keys;
  int isortheader;
  int sizeof_my_indxkey;
  bool allzero;
  bool passrangetest;
  */
  char* ilinekey;
  char* xlinekey;
  int ilineinc;
  int indx_ilinekey;
  int indx_xlinekey;
  float ilinemin;
  float ilinemax;
  float xlinemin;
  float xlinemax;
  int headers_in_buffer;
  int i_header;
  int ii_header;
  int n_headers_read;
  int n_headers_left;
  double doubleiline;
  float** header_buf=NULL;
  int indx_number;
  int min_indx_number;
  int max_indx_number;

  
  sf_init (argc,argv);

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls amount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
  /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  if(verbose>0)fprintf(stderr,"verbose=%d\n",verbose);
 
  /*****************************************/
  /* initialize the input and output files */
  /*****************************************/
  if(verbose>0)fprintf(stderr,"read name of input file name\n");
  
  infile_filename=sf_getstring("input");
  /* \n
     Input file for traces amplitudes
  */
  if(infile_filename==NULL) infile = sf_input ("in");
  else infile = sf_input (infile_filename);

  if(verbose>0)
    fprintf(stderr,"set up output file for tah - should be stdout\n");
  out = sf_output ("out");

  if(verbose>0)fprintf(stderr,"read name of input headers file\n");
  headers_filename=sf_getstring("headers");
  /* \n
     Trace header file name.  Default is the input data file
     name, with the final .rsf changed to _hdr.rsf 
  */

  if(headers_filename==NULL){
    /* compute headers_filename from infile_filename by replacing the final
       .rsf with _hdr.rsf */
    if(!(0==strcmp(infile_filename+strlen(infile_filename)-4,".rsf"))){
	fprintf(stderr,"parameter input, the name of the input file, does\n");
	fprintf(stderr,"not end with .rsf, so header filename cannot be\n");
	fprintf(stderr,"computed by replacing the final .rsf with _hdr.rsf.\n");
	sf_error("default for headers parameter cannot be computed.");
    }
    /* allocate header file name to be length of infile_filename +4+1 to allow 
       _hdr to be inserted and \0 to be appended.  */
    headers_filename=malloc(strlen(infile_filename)+4+1);
    headers_filename[0]='\0'; /* make string zero length, null terminated */
    strcpy(headers_filename,infile_filename);
    strcpy(headers_filename+strlen(infile_filename)-4,"_hdr.rsf\0");
    if(verbose>2)
      fprintf(stderr,"parameter header defaulted.  Computed to be #%s#\n",
			 headers_filename);
  }
  if(verbose>2)
    fprintf(stderr,"parameter header input or computed  #%s#\n",
		       headers_filename);

  inheaders = sf_input(headers_filename);

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
  if(verbose>1)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  header = sf_floatalloc(n1_headers);
 
  /* I make stdout just trace header and one trace amplitude for now.
     I do not think it will have a use and will delete it later */
  if(verbose>1)
    fprintf(stderr,"need to add 2 words, record type and record length.\n");
  sf_putint(out,"n1",1+n1_headers+2);
  sf_putint(out,"n1_traces",1);

  /*************************************************************************/
  /* put info into history so later tah programs can figure out the headers*/
  /*************************************************************************/
  if(verbose>1)fprintf(stderr,"put n1_headers to history\n");
  sf_putint(out,"n1_headers",n1_headers);
  if(verbose>1)fprintf(stderr,"test typehead\n");
  if(SF_INT == typehead) sf_putstring(out,"header_format","native_int");
    else                 sf_putstring(out,"header_format","native_float");
  sf_putstring(out,"headers",headers_filename);
  /* the list of any extra header keys */
  if(verbose>1)fprintf(stderr,"segy_init\n");
  segy_init(n1_headers,inheaders);
  if(verbose>1)fprintf(stderr,"segy2hist\n");
  segy2hist(out,n1_headers);
  
  /* put the history from the input file to the output */
  if(verbose>1)fprintf(stderr,"fileflush out\n");
  sf_fileflush(out,infile);
  
  /* get the index keys.  Normally this will be iline and xline */
  if(!(ilinekey=sf_getstring("iline")))ilinekey="iline";
  /* \n
     header key for the main index key.  This should be iline, but you 
     may have non-standard trace headers or a wierd use of this program 
  */
  if(!(xlinekey=sf_getstring("xline")))xlinekey="xline";
  /* \n
     header key for the secondary index key.  This should be xline, but you
     may have non-standard trace headers or a wierd use of this program 
  */
  
  /* ilinemin, ilinemax, xlinemin, and xlinemaxcan be used to remove
     be trace headers with nulls. They default to -1e31, +1e31 */
  if(!sf_getfloat("ilinemin",&ilinemin))ilinemin=-1e31;
  /* \n
     minimum "iline" header key to include in the index.  Use this parameter
     to remove null trace headers or traces outside project area.
  */
  if(!sf_getfloat("ilinemax",&ilinemax))ilinemax=-1e31;
  /* \n
     maximum "iline" header key to include in the index.  Use this parameter
     to remove null trace headers or traces outside project area.
  */
  if(!sf_getfloat("xlinemin",&xlinemin))xlinemin=-1e31;
  /* \n
     minimum "xline" header key to include in the index.  Use this parameter
     to remove null trace headers or traces outside project area.
  */
  if(!sf_getfloat("xlinemax",&xlinemax))xlinemax=-1e31;
  /* \n
     maximum "xline" header key to include in the index.  Use this parameter
     to remove null trace headers or traces outside project area.
  */

  sf_getstring("indxdir");
  /* \n
     The name of the directory containing the iline,xline index.  This 
     directory will be in DATAPATH (probably the environment variable). The 
     directory also continues a file "filenames", a list of the trace and 
     header files that contributes to this index. The directory contains files 
     with names "indx#" here # is an integer multiple of ilineinc. These files 
     contains a record for each contributing trace with filenumber, 
     tracenumber, and the trace header. The file containing the trace is 
     determined using the  can be read by using the filenumber and the 
     "filenames" file.  The tracenumber defines the location of the trace 
     in the file.
  */

  if(!sf_getint("ilineinc",&ilineinc))ilineinc=10;
  /* \n
     incrment in iline for the index
  */

  indx_ilinekey=segykey(ilinekey);
  indx_xlinekey=segykey(xlinekey);

  fprintf(stderr,"process the filenames file\n");  
  /* kls process the file, "filenames"
     -if it does not exist create it
     -if it already exists, count the number files already in "filenemae"
     -add the traceheaders and traces files to "filenames"
     -remenber the file_number of this trace/header input file
  */

  /* make buffer 1 gbyte, about 25 million traces. that should be reasonable */
  headers_in_buffer=1000*1000*1000/(4*n1_headers);
  header_buf=sf_floatalloc2(n1_headers,headers_in_buffer);
  /* sf_alloc2 for array of trace headers */
  if(verbose>0)fprintf(stderr,"start loop to read headers to build sort key\n");
  for (i_header=0; i_header<n_headers; i_header+=headers_in_buffer){
    /* read smaller of headers_in_buffer and number of headers left to read */
    fprintf(stderr,"read a buffer of trace headers\n");
    n_headers_read=headers_in_buffer;
    n_headers_left=n_headers-i_header;
    if(n_headers_read>n_headers_left)n_headers_read=n_headers_left;
    fprintf(stderr,"read words=%d\n",n1_headers*n_headers_read);
    sf_floatread(&header_buf[0][0],n1_headers*n_headers_read,inheaders);
    fprintf(stderr,"read the headers\n");
    /* find the minimum and maximum indx_number this buffer contributes to */
    if(typehead==SF_FLOAT)doubleiline=        header_buf [0][indx_ilinekey];
    else                  doubleiline=((int**)header_buf)[0][indx_ilinekey];
    indx_number=rounddouble(doubleiline/((double)ilineinc))*ilineinc;
    min_indx_number=max_indx_number=indx_number;
    fprintf(stderr,"loop to find min/max\n");
    for (ii_header=0; ii_header<n_headers_read; ii_header++){
      if(typehead==SF_FLOAT){
	doubleiline=        header_buf [ii_header][indx_ilinekey];
      }else{                  
	doubleiline=((int**)header_buf)[ii_header][indx_ilinekey];
      }
      fprintf(stderr,"ii_header=%d,doubleiline=%f %d\n",
              ii_header,doubleiline,
	      ((int**)header_buf)[ii_header][indx_ilinekey]);
      indx_number=rounddouble(doubleiline/((double)ilineinc))*ilineinc;
      if(min_indx_number>indx_number)min_indx_number=indx_number;
      if(max_indx_number<indx_number)max_indx_number=indx_number;
    }
    fprintf(stderr,"min_indx_number=%d, max_indx_number=%d\n",
	    min_indx_number   ,max_indx_number);
    /* kls pass over each indx_file:
       -open indx_file (create it if it does not exist)
       -loop over all the headers in header_buf
          -write filenumber, traceindx, and the trace header to the end of
	     indx_file 
       -close the indx_file   */
    /* to convert relative to absolute path:
       #include <stdio.h>
       #include <stdlib.h>
       int main(int c,char**v) {
  
         fprintf(stderr,"%s\n",realpath(v[1],0));
         return 0;
       }

    */







  } /* end of the for loop that reads a buffer of headers */









  exit(0);
}

  
