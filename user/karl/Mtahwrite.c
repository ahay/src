/* Read Trace And Header (tah) from standard input, write to separate files

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

The sftahwrite program reads the trace and header data (tah) from 
standard input (usually a pipe), separates the trace data from the
header data.  The trace data is written to output and the header is
written to outheaders.  The trace headers are used to select the 
location in the file to write the data.  The iline and xline headers 
are used in the following example to put stacked data in 
(time, xline, iline) order so it can be viewed using sfgrey.

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
   10/22/2013 Karl Schleicher Factor reused functions into tahsub.c
   01/10/2014 Karl Schleicher Provide headers default, documentation 
*/
#include <string.h>
#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

/* very sparingly make some global variables. */
int verbose;

int main(int argc, char* argv[])
{
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  int dim;
  off_t n_in[SF_MAX_DIM];
  int iaxis;
  int dim_output;
  int *indx_of_keys;
  bool label_argparmread,n_argparmread,o_argparmread,d_argparmread;
  char parameter[13];
  char* label[SF_MAX_DIM];
  off_t n_output[SF_MAX_DIM];
  off_t n_outheaders[SF_MAX_DIM];
  float o_output[SF_MAX_DIM];
  float d_output[SF_MAX_DIM];
  sf_axis output_axa_array[SF_MAX_DIM];
  sf_file output=NULL, outheaders=NULL;
  char* output_filename=NULL;
  char* outheaders_filename=NULL;
  
  sf_init (argc,argv);

   /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
  /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  fprintf(stderr,"verbose=%d\n",verbose);
 
  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/

  if(verbose>0)fprintf(stderr,"read infile name\n");  
  in = sf_input ("in");

  if(verbose>0)
    fprintf(stderr,"read outfile name (probably default to stdout\n");
  out = sf_output ("out");

  if (!sf_histint(in,"n1_traces",&n1_traces))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n1_headers",&n1_headers)) 
    sf_error("input data does not define n1_headers");

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  fheader = sf_floatalloc(n1_headers);

  fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  
/* maybe I should add some validation that n1== n1_traces+n1_headers+2
     and the record length read in the second word is consistent with 
     n1.  */

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************************/
  /* set up the output and outheaders files for the traces and headers */
  /********************************************************************/

  output_filename=sf_getstring("output");
  /* \n
     output trace filename. Required parameter with no default.
  */
  if(NULL==output_filename) sf_error("output is a required parameter");
  output=sf_output(output_filename);

  outheaders_filename=sf_getstring("outheaders");
  /* \n
     Output trace header file name.  Default is the input data
     file name, with the final .rsf changed to _hdr.rsf. 
  */  

  if(outheaders_filename==NULL){
    /* compute headers_filename from output_filename by replacing the final
       .rsf with _hdr.rsf */
    if(!(0==strcmp(output_filename+strlen(output_filename)-4,".rsf"))){
	fprintf(stderr,"parameter output, the name of the output file,\n");
	fprintf(stderr,"does not end with .rsf, so header filename cannot\n");
	fprintf(stderr,"be computed by replacing the final .rsf with\n");
	fprintf(stderr,"_hdr.rsf.\n");
	sf_error("default for outheaders parameter cannot be computed.");
    }
    outheaders_filename=malloc(strlen(output_filename)+60);
    strcpy(outheaders_filename,output_filename);
    strcpy(outheaders_filename+strlen(output_filename)-4,"_hdr.rsf\0");
    if(verbose>1)
      fprintf(stderr,"parameter outheader defaulted.  Computed to be #%s#\n",
			 outheaders_filename);
  }
  if(verbose>1)fprintf(stderr,"parameter outheader input or computed  #%s#\n",
		       outheaders_filename);
  outheaders=sf_output(outheaders_filename);
 
  /* get each of the axis information: 
     label2, n2, o2, d2,  
     label3, n3, o3, d3,
     etc
     label1, n1, o1, d1  is always defaulted from input */
  sf_putint   (output    ,"n1"    ,n1_traces );
  sf_putint   (outheaders,"n1"    ,n1_headers);
  sf_putfloat (outheaders,"d1"    ,1         );
  sf_putfloat (outheaders,"o1"    ,0         );
  sf_putstring(outheaders,"label1","none"    );
  sf_putstring(outheaders,"unit1" ,"none"    );
  
  dim_output=1;
  for (iaxis=1; iaxis<SF_MAX_DIM; iaxis++){
    label_argparmread=n_argparmread=o_argparmread=d_argparmread=false;
    sprintf(parameter,"label%d",iaxis+1);
    fprintf(stderr,"try to read %s\n",parameter);
    if ((label[iaxis]=sf_getstring(parameter))) {
      /*(label#=(2,...)  name of each of the axes. 
	   label1 is not changed from input. Each label must be a 
	   header key like cdp, cdpt, or ep.  The trace header 
	   values are used to define the output trace location in
	   the output file. )*/
      fprintf(stderr,"got %s=%s\n",parameter,label[iaxis]);
      sf_putstring(output    ,parameter,label[iaxis]);
      sf_putstring(outheaders,parameter,label[iaxis]);
      label_argparmread=true;
    }
    sprintf(parameter,"n%d",iaxis+1);
    fprintf(stderr,"try to read %s\n",parameter);
    if (sf_getlargeint  (parameter,&n_output[iaxis])) {
      /*( n#=(2,...) number of locations in the #-th dimension )*/ 
      fprintf(stderr,"got %s=%lld\n",parameter,(long long) n_output[iaxis]);
      sf_putint(output    ,parameter,n_output[iaxis]);
      sf_putint(outheaders,parameter,n_output[iaxis]);
      n_argparmread=true;
    }
    sprintf(parameter,"o%d",iaxis+1);
    if (sf_getfloat(parameter,&o_output[iaxis])) {
      /*( o#=(2,...) origin of the #-th dimension )*/ 
      sf_putfloat(output    ,parameter,o_output[iaxis]);
      sf_putfloat(outheaders,parameter,o_output[iaxis]);
      o_argparmread=true;
    }
    sprintf(parameter,"d%d",iaxis+1);
    if (sf_getfloat(parameter,&d_output[iaxis])) {
      /*( d#=(2,...) delta in the #-th dimension )*/ 
      sf_putfloat(output    ,parameter,d_output[iaxis]);
      sf_putfloat(outheaders,parameter,d_output[iaxis]);
      d_argparmread=true;
    }
    if(!label_argparmread && !n_argparmread && 
       !    o_argparmread &&  !d_argparmread){
      /* none of the parameter were read
	 you read all the parameters in the previous iteration
	 compute the output dimension and exit the loop */
      dim_output=iaxis;
      break;
    }
    if(label_argparmread && n_argparmread && o_argparmread && d_argparmread){
      /* all the parameters for thisi axis were read.  loop for next iaxis */
      if(verbose>0){
	fprintf(stderr,"label, n, i, and d read for iaxis%d\n",iaxis+1);
      }
    } else {
      sf_warning("working on iaxis=%d. Program expects to read\n",iaxis+1);
      sf_warning("label%d, n%d, i%d, and d%d from command line.\n",
		 iaxis+1,iaxis+1,iaxis+1,iaxis+1);
      if(!label_argparmread) sf_error("unable to read label%d\n",iaxis+1); 
      if(!n_argparmread    ) sf_error("unable to read n%d    \n",iaxis+1); 
      if(!o_argparmread    ) sf_error("unable to read o%d    \n",iaxis+1); 
      if(!d_argparmread    ) sf_error("unable to read d$d    \n",iaxis+1); 
    }
  }
  
  /* if the input file is higher dimention than the output, then the 
     size of all the higher dimensions must be set to 1. */
  /* kls use sf_axis to get structure for each axis with n, o, d, l, and u
     this can be used to compute the index for sf_seek */
  dim = sf_largefiledims(in,n_in);
  for (iaxis=dim_output; iaxis<dim; iaxis++){
    sprintf(parameter,"n%d",iaxis+1);
    sf_putint(output,parameter,n_output[iaxis]);
    sf_putint(outheaders,parameter,n_output[iaxis]);
  }
  /* for a test zero n_output and dim_output and see it you can read the 
     history file */

  sf_largefiledims(output,n_output);
  sf_largefiledims(outheaders,n_outheaders);
  if(verbose>1){
    for (iaxis=0; iaxis<SF_MAX_DIM; iaxis++){
      fprintf(stderr,"from sf_largefiledims(output.. n%d=%lld outheaders=%lld\n",
	      iaxis+1,(long long) n_output[iaxis],(long long) n_outheaders[iaxis]);
    }
  }

  
  /* the list header keys (including any extras) and get the index into 
     the header vector */
  segy_init(n1_headers,in);
  indx_of_keys=sf_intalloc(dim_output);
  for (iaxis=1; iaxis<dim_output; iaxis++){
    /* kls need to check each of these key names are in the segy header and
       make error message for invalid keys.  Of does segykey do this? NO, just
       segmentation fault. */
    if(verbose>0)fprintf(stderr,"get index of key for %s\n",label[iaxis]); 
    indx_of_keys[iaxis]=segykey(label[iaxis]);
  }

  sf_fileflush(output,in);

  sf_putint(outheaders,"n1",n1_headers);
  if(0==strcmp("native_int",sf_histstring(in,"header_format"))){
    sf_settype(outheaders,SF_INT);
  } else {
    sf_settype(outheaders,SF_FLOAT);
  }
  sf_fileflush(outheaders,in);

  fprintf(stderr,"start trace loop\n");

  /* kls maybe this should be in function sf_tahwritemapped_init */

  {
    off_t file_offset;
    off_t i_output[SF_MAX_DIM];
    float temp_float=0.0;

    for(iaxis=0; iaxis<SF_MAX_DIM; iaxis++){
      i_output[iaxis]=n_output[iaxis]-1;
    }
    file_offset=sf_large_cart2line(dim_output,n_output,i_output)*sizeof(float);
    sf_seek(output,file_offset,SEEK_SET);
    sf_floatwrite(&temp_float,1,output);
    i_output[0]=n_outheaders[0]-1;
    file_offset=sf_large_cart2line(dim_output,n_outheaders,i_output)*
                      sizeof(float);
    sf_seek(outheaders,file_offset,SEEK_SET);
    sf_floatwrite(&temp_float,1,outheaders);
  }

  for (iaxis=0; iaxis<SF_MAX_DIM; iaxis++){
      /* sf_axis temp; */
    output_axa_array[iaxis]=sf_iaxa(output,iaxis+1);
    if(verbose>2){
	/* temp=output_axa_array[iaxis]; */
      fprintf(stderr,"axis=%d sf_n(output_axa_array[iaxis])=%d\n",
	              iaxis+1,sf_n(output_axa_array[iaxis]));
    }
    /* kls why does this fail? 
       fprintf(stderr,"temp->n=%d\n",temp->n);
    */
  }

  /***************************/
  /* start trace loop        */
  /***************************/
  fprintf(stderr,"start trace loop\n");
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1){
      for(iaxis=2; iaxis<dim_output; iaxis++){
	fprintf(stderr,"label[%d]=%s",
		iaxis,label[iaxis]);
	if(typehead == SF_INT)fprintf(stderr,"%d",
				      ((int*)fheader)[indx_of_keys[iaxis]]);
	else                  fprintf(stderr,"%f",
				      fheader[indx_of_keys[iaxis]]);
      }
      fprintf(stderr,"\n");
    }
    tahwritemapped(verbose,intrace, fheader, 
		   n1_traces, n1_headers,
		   output, outheaders,
		   typehead, output_axa_array,
		   indx_of_keys, dim_output,
		   n_output,n_outheaders);
    /**********************************************/
    /* write trace and headers to the output pipe */
    /**********************************************/
    put_tah(intrace, fheader, n1_traces, n1_headers, out);
  }

  exit(0);
}

  
