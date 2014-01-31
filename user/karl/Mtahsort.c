/* Read Trace And Header from separate files in sorted order, write to pipe

tah is the abbreviation of Trace And Header.  Madagascar programs 
that begin with sftah are designed to:
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
sftahwrite writes the trace data to mappedstack.rsf and the headers 
are written to the file mappedstack_hdr.rsf.  The order of the data in
the output file is defined by the iline and xline trace headers, so the 
data order is (time,xline,iline).  Finally, the output volume is
displayed using sfgrey.
*/
/*
   Program change history:
   date       Who             What
   01/19/2014 Karl Schleicher Original program based in Mtahread
*/
#include <string.h>

#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

/* sparingly make some global variables. */
int verbose;
int numsortheadernames;

typedef struct {
  int filenumber;
  int tracenumber;
  double keyvalue[1];
} sortkey;

int qsort_compare(const void *a, const void *b){
  int ikey;
  sortkey* a_sortkey=*(sortkey**)a;
  sortkey* b_sortkey=*(sortkey**)b;

  for (ikey=0; ikey<numsortheadernames; ikey++){
    /* as soon as keys differ you can return the compare */
    if(a_sortkey->keyvalue[ikey] > b_sortkey->keyvalue[ikey])return 1;
    if(a_sortkey->keyvalue[ikey] < b_sortkey->keyvalue[ikey])return -1;
  }
  /* keys are all equal.  */
  return 0;
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
  float* intrace=NULL;
  int i_trace;
  char* infile_filename=NULL;
  char* headers_filename=NULL;
  sortkey** myarray_pointers_to_sortkey;
  char* sort;
  char** sortheadernames;
  char** sorttokens1;
  double* signsortname;
  double* sortmin;
  double* sortmax;
  double* sortinc;
  int numtracestosort;
  int size_myarray_pointers_to_sortkey;
  int i1_headers;
  int* indx_of_keys;
  int isortheader;
  int bytes_in_my_sortkey;
  bool allzero;
  bool passrangetest;
  
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
    headers_filename=malloc(strlen(infile_filename)+60);
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
 
  if(verbose>1)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);

  intrace= sf_floatalloc(n1_traces);

  if(verbose>1)
    fprintf(stderr,"need to add 2 words, record type and record length.\n");
  sf_putint(out,"n1",n1_traces+n1_headers+2);
  sf_putint(out,"n1_traces",n1_traces);

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

  /* get the sort key */
  sort=sf_getstring("sort");
  /* /n
     list of the sort keys.  Each key must be a trace header key name.
     It may be preceeded with + (the default) for ascending or - for 
     descending sort direction.  The key may be followed with :min,max 
     or :min,max,inc.  These parameters allow you to select a subset of 
     of the traces based on the header key.  The sort keys are blank
     seperated you should enclose the sort string in \"\'s.  Examples are
     sort=\"iline xline offset\" and sort=\"cdp:100,500,25 offset\"
  */
  if(NULL==sort)sf_error("sort not input\n");
  sorttokens1=delimittedtokens2list(sort," ",&numsortheadernames);
  if(verbose>0){
    for (isortheader=0; isortheader<numsortheadernames; isortheader++){
      fprintf(stderr,"sorttokens1[%d]=%s\n",
	              isortheader,sorttokens1[isortheader]);
    }
  }

  sortheadernames=(char**)sf_alloc(numsortheadernames,sizeof(char*)); 
  signsortname=   (double*)sf_alloc(numsortheadernames,sizeof(double));
  sortmin     =   (double*)sf_alloc(numsortheadernames,sizeof(double));
  sortmax     =   (double*)sf_alloc(numsortheadernames,sizeof(double));
  sortinc     =   (double*)sf_alloc(numsortheadernames,sizeof(double));

  for (isortheader=0; isortheader<numsortheadernames; isortheader++){
    char** sorttokenparts;
    int numsorttokenparts;
    char** pointer_to_minmaxinc;
    int num_minmaxinc;
    /* look for +/- at the beginning of the token and use it to set
       signsortname[isortheader].  This is used to sort ascending or
       desceinding */
    if('-'==sorttokens1[isortheader][0])signsortname[isortheader]=-1.0;
    else                                signsortname[isortheader]=1.0;
    if('+'==sorttokens1[isortheader][0] || '-'==sorttokens1[isortheader][0]){
      /* first char is + or -.  Remove it by increment the pointer by 1. */
      (sorttokens1[isortheader])++;
    }
    fprintf(stderr,"trimmed sorttokens1[%d]=%s\n",
	        isortheader,sorttokens1[isortheader]);
    /* the token now may be a header:min,max,inc (eg cdp:100,500,25).  look 
       for a :.  If there is one, use it to split the header name from the 
       min,max,inc */
    sorttokenparts=delimittedtokens2list(sorttokens1[isortheader],":",
					 &numsorttokenparts);
    /* default min,max,inc.  Then look for input values */
    sortmin[isortheader]=-1e307; /* Is range of doubles +- 1.7e308 ? */ 
    sortmax[isortheader]= 1e307; /* Is range of doubles +- 1.7e308 ? */ 
    sortinc[isortheader]=-1;     /* negative inc means not input */
    if(numsorttokenparts>1){
      /* there is a : in this sort key.  If there are more than one :
	 there is a syntax error.  Report and exit */
      if(numsorttokenparts>2){
	sf_error("error reading sort key %s.  There is more than one :.\n",
		 sorttokens1[isortheader]);
      }
      fprintf(stderr,"sorttokenpart[1]=%s\n",sorttokenparts[1]);
      /* now dig the min,max,inc out kls */
      pointer_to_minmaxinc=delimittedtokens2list(sorttokenparts[1],",",
						 &num_minmaxinc);
      if(num_minmaxinc !=2 && num_minmaxinc !=3){
	sf_warning("parsing minmaxinc=%s.\n",sorttokenparts[1]);
	sf_error("This should be 2 or 3 floats seperated by commas\n");
      }
      /* kls is this bullet proof? */
      sortmin[isortheader]=strtod(pointer_to_minmaxinc[0],NULL);
      sortmax[isortheader]=strtod(pointer_to_minmaxinc[1],NULL);
      if(num_minmaxinc ==3)
	sortinc[isortheader]=strtod(pointer_to_minmaxinc[2],NULL);
    }	 
    /* save sorttokenparts[0],  the portion of sorttokens1[isortheader] 
       without initial +/- and without the stuff after : */
    sortheadernames[isortheader]=sorttokenparts[0];
    if(verbose>1)
      fprintf(stderr,"copied sortheadernames[%d]=%s min=%e max=%e inc=%e\n",
                      isortheader,sortheadernames[isortheader],
                      sortmin[isortheader],sortmax[isortheader],
	              sortinc[isortheader]);
    /* need to hang onto sorttokenparts, sortheadernames[isortheader] 
       points to a string in that structure. do not free.  Let little
       bits of structure memory leak. */
  }
  if(verbose>0){
    for (isortheader=0; isortheader<numsortheadernames; isortheader++){
      fprintf(stderr,"sortheadernames[%d]=%s min=%e max=%e inc=%e\n",
	              isortheader,sortheadernames[isortheader],
	              sortmin[isortheader],sortmax[isortheader],
                      sortinc[isortheader]);
    }
  }
  indx_of_keys=sf_intalloc(numsortheadernames);
  for(isortheader=0; isortheader<numsortheadernames; isortheader++){
    indx_of_keys[isortheader]=segykey(sortheadernames[isortheader]);
  }
  numtracestosort=0;
  size_myarray_pointers_to_sortkey=100;
  bytes_in_my_sortkey=2*sizeof(*myarray_pointers_to_sortkey)+
                        sizeof(double)*(numsortheadernames-1);
  myarray_pointers_to_sortkey=malloc(size_myarray_pointers_to_sortkey*
                                     sizeof(myarray_pointers_to_sortkey));

  if(verbose>0)fprintf(stderr,"start loop to read headers to build sort key\n");
  for (i_trace=0; i_trace<n_traces; i_trace++){
    if(verbose>2 ||(verbose>0 && i_trace<5)){
      fprintf(stderr,"i_trace=%d\n",i_trace);
    }
    /****************/
    /* read headers */
    /****************/
    sf_floatread(header,n1_headers,inheaders);
    /* if trace header is all zero, skip trace using break from loop */
    allzero=true;
    for(i1_headers=0; i1_headers<n1_headers; i1_headers++){
      if(0!=header[i1_headers]){
	allzero=false;
	break;
      }
    }
    if(allzero)continue; /* do not add header to list of keys */

    /* test if trace passes the sort ranges */
    passrangetest=true;
    for(isortheader=0; isortheader<numsortheadernames; isortheader++){
      double myvalue;
      if(typehead==SF_FLOAT)
	myvalue=header[indx_of_keys[isortheader]];
      else
	myvalue=((int*)header)[indx_of_keys[isortheader]];
      /* add test for inc */
      if(myvalue<sortmin[isortheader] || 
	 myvalue>sortmax[isortheader]){
	passrangetest=false;
	break;
      }
      if(sortinc[isortheader]>0){
	double delta=(myvalue-sortmin[isortheader])/sortinc[isortheader];
	if(1e-5<fabs(delta-llround(delta))){
	  passrangetest=false;
	  break;
	}
      }
    }
    if(false==passrangetest) continue; /* do not add header to list of keys
                                          header is outsidd the sort ranges */
			       
    /* add this trace to the headers to sort */
    if(numtracestosort>=size_myarray_pointers_to_sortkey){
      /* reallocate the array */
      size_myarray_pointers_to_sortkey*=1.2;
      myarray_pointers_to_sortkey=realloc(myarray_pointers_to_sortkey,
	      size_myarray_pointers_to_sortkey*
	         sizeof(myarray_pointers_to_sortkey));
    }
    myarray_pointers_to_sortkey[numtracestosort]=malloc(bytes_in_my_sortkey);
    for(isortheader=0; isortheader<numsortheadernames; isortheader++){
      if(typehead==SF_FLOAT)
	myarray_pointers_to_sortkey[numtracestosort]->keyvalue[isortheader]=
	header[indx_of_keys[isortheader]];
      else
	myarray_pointers_to_sortkey[numtracestosort]->keyvalue[isortheader]=
	  ((int*)header)[indx_of_keys[isortheader]];
      /* multiply header value by signsortkey to handle sort in decreasing 
	 order */
      myarray_pointers_to_sortkey[numtracestosort]->keyvalue[isortheader]*=
	signsortname[isortheader];
    }
    myarray_pointers_to_sortkey[numtracestosort]->tracenumber=i_trace;
    myarray_pointers_to_sortkey[numtracestosort]->filenumber=0;

    numtracestosort++;
  }
  qsort(myarray_pointers_to_sortkey,numtracestosort,
	sizeof(myarray_pointers_to_sortkey),qsort_compare);
  /* print the headers */
  if(verbose>3){
    for(i_trace=0; i_trace<numtracestosort; i_trace++){
      for(isortheader=0; isortheader<numsortheadernames; isortheader++){
	fprintf(stderr,"%s = %f ",sortheadernames[isortheader],
		myarray_pointers_to_sortkey[i_trace]->keyvalue[isortheader]);
      }
      fprintf(stderr,"tracenumber=%d ",
	      myarray_pointers_to_sortkey[i_trace]->tracenumber);
      fprintf(stderr,"filenumber=%d ",
	      myarray_pointers_to_sortkey[i_trace]->filenumber);
      fprintf(stderr,"\n");
    }
  }

  if(verbose>0)fprintf(stderr,"start read trace loop n_traces=%d\n",n_traces);
  for (i_trace=0; i_trace<numtracestosort; i_trace++){
    off_t file_tracenum;
    if(verbose>2 ||(verbose>0 && i_trace<5)){
      fprintf(stderr,"i_trace=%d\n",i_trace);
    }

    /**************************/
    /* read trace and headers */
    /**************************/
    file_tracenum=myarray_pointers_to_sortkey[i_trace]->tracenumber;

    sf_seek(inheaders,  file_tracenum*n1_headers*sizeof(float),  SEEK_SET);
    sf_floatread(header,n1_headers,inheaders);

    sf_seek(infile,  file_tracenum*n1_traces*sizeof(float),  SEEK_SET);
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

  
