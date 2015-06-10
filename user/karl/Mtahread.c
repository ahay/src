/* Read Trace And Header from separate files, combine, write to pipe

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
   04/26/2013 Karl Schleicher Original program
   01/10/2014 Karl Schleicher Provide headers default, documentation
   06/05/2015 Karl Schleicher add makeheader=y to load axis coords to hdrs
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
  bool makeheader;
  int iaxis;
  char parameter[13];
  char* label_in[SF_MAX_DIM+1];
  int n_in[SF_MAX_DIM+1];
  float o_in[SF_MAX_DIM+1];
  float d_in[SF_MAX_DIM+1];
  int indx_of_keys[SF_MAX_DIM+1];
  int traceindx[SF_MAX_DIM+1];
  float tracecoord[SF_MAX_DIM+1];
  int indxleft;
  
  sf_init (argc,argv);

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
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
     Input file for traces amplitudes.  You can list a file here, has the 
     input file name will be used to compute the name of the header files.

     The input trace amplitudes can also be read from standard input by
     just supplying standard input and omitting this paramater,  This 
     is useful it you wish to do sopme initial processing of the input
     rsf file containing the trace amplitudes.  This is useful if you need 
     to change input axis labels to use the makeheader=yes.
  */
  if(infile_filename==NULL) infile = sf_input ("in");
  else infile = sf_input (infile_filename);

  /* get the axis labels on the input file */
  /* these will be used if users requests makeheaders */
  for (iaxis=2; iaxis<SF_MAX_DIM+1; iaxis++){
    sprintf(parameter,"label%d",iaxis);
    if(verbose>1)fprintf(stderr,"first try to read %s\n",parameter);
    /* get label, o, n, and i for each axis.  These will be used to 
       make headers from axis (if that option is requested) */
    if (!(label_in[iaxis]=sf_histstring(infile,parameter))) {
      label_in[iaxis]="none";
    }
    sprintf(parameter,"o%d",iaxis);
    if(verbose>1)fprintf(stderr,"first try to read %s\n",parameter);
    if (!sf_histfloat(infile,parameter,&o_in[iaxis])){
      label_in[iaxis]="none";
    }
    sprintf(parameter,"d%d",iaxis);
    if(verbose>1)fprintf(stderr,"first try to read %s\n",parameter);
    if (!sf_histfloat(infile,parameter,&d_in[iaxis])){
      label_in[iaxis]="none";
    }
    sprintf(parameter,"n%d",iaxis);
    if(verbose>1)fprintf(stderr,"first try to read %s\n",parameter);
    if (!sf_histint(infile,parameter,&n_in[iaxis])){
      n_in[iaxis]=1;
    }
    if(verbose>1){
      fprintf(stderr,"on axis=%d got label,o,d,n=%s,%f,%f,%d \n",iaxis,
	      label_in[iaxis],o_in[iaxis],d_in[iaxis],n_in[iaxis]);
    }
  }
  /* compute index_of_keys after calling segy_init */ 






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
    if(verbose>2){
      fprintf(stderr,"parameter header defaulted.  Computed to be #%s#\n",
	      headers_filename);
    }
  }
  if(verbose>2)fprintf(stderr,"parameter header input or computed  #%s#\n",
		       headers_filename);

  if(!sf_getbool("makeheader",&makeheader))makeheader=false;
  /* \n
     Option to load headers using the input file axis labels.  If axis 
     label2 through label9 match a header key then that coordinate is
     loaded to the traces header.  This can be used to load the source
     coordinate to the sx header location.  This may require changing
     the axis label because Madagascar axis labels are not the same as
     segy trace headers.  for example axis 2 coordiante can be loaded in
     sx trace header by:
        <spike.rsf sfput label2=sx \\
           | sftahread headers=spike_hdr.rsf makeheader=y \\ 
           | sftahgethw key=sx >/dev/null
 */
  
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

  if(verbose>1)fprintf(stderr,"need to add 2 words.  Record type and record length\n");
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
  if(0)  {
    int tempint;
    tempint=segykey("ep");
    fprintf(stderr,"ep tempint=%d\n",tempint);
    tempint=segykey("epx");
    fprintf(stderr,"epx tempint=%d\n",tempint);
  }
  if(makeheader){
    for (iaxis=2; iaxis<SF_MAX_DIM+1; iaxis++){
      if(0==strcmp("none",label_in[iaxis])){
	indx_of_keys[iaxis]=-1;
      } else {
	indx_of_keys[iaxis]=segykey(label_in[iaxis]);
	if(indx_of_keys[iaxis]<0){
	  sf_warning("************************************************");
	  sf_warning("************************************************");
	  sf_warning("axis %d has label that is not a header key",iaxis);
	  sf_warning("This axis label will not be loaded");
	  sf_warning("************************************************");
	  sf_warning("************************************************");
	}
      }
      if(verbose>1){
	fprintf(stderr,"indx_of_keys[%d]=%d\n",iaxis,indx_of_keys[iaxis]);
      }
    }
  }
  /* put the history from the input file to the output */
  if(verbose>1)fprintf(stderr,"fileflush out\n");
  sf_fileflush(out,infile);

  if(verbose>0)fprintf(stderr,"start trace loop n_traces=%d\n",n_traces);
  for (i_trace=0; i_trace<n_traces; i_trace++){
    if(verbose>2 ||(verbose>0 && i_trace<5)){
      fprintf(stderr,"i_trace=%d\n",i_trace);
    }
    /**************************/
    /* read trace and headers */
    /**************************/
    sf_floatread(header,n1_headers,inheaders);
    sf_floatread(intrace,n1_traces,infile);

    if(makeheader){
      indxleft=i_trace;
      for (iaxis=2; iaxis<SF_MAX_DIM+1; iaxis++){
	traceindx[iaxis]=indxleft%n_in[iaxis];
	indxleft/=n_in[iaxis];
	tracecoord[iaxis]=o_in[iaxis]+traceindx[iaxis]*d_in[iaxis];
	if(verbose>1){
	  fprintf(stderr,"i_trace=%d,iaxis=%d,tracecoord[iaxis]=%f\n",
		  i_trace,   iaxis,   tracecoord[iaxis]);
	}
	if(indx_of_keys[iaxis]>=0){
	  /* kls what's needed to add to make this work with integer headers?*/
	  header[indx_of_keys[iaxis]]=tracecoord[iaxis];
	}
      }
    }
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

  
