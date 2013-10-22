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
#include <rsfsegy.h>

char** sf_getstringarray(char* key, int* numkeys){

  char** list_of_keys;
  char* par;
  char* one_par;
  char* copy_of_par;
  int ikey;

  par=sf_getstring(key);
  if(par==NULL){
    return NULL;
  }

  /* count number of keys */
  *numkeys=0;
  copy_of_par=sf_charalloc(strlen(par));
  strcpy (copy_of_par, par);
  one_par=strtok(copy_of_par,",");
  while (one_par!=NULL){
    (*numkeys)++;
    one_par=strtok(NULL,",");
  }
  /* break the long string into an array  of strings */
  list_of_keys=(char**)sf_alloc(*numkeys,sizeof(char*));
  /* no need to reallocate, but including free(copy_of_par) clobberred part of
    list_of_keys.  This spooked me on exactly how strtok works */
  copy_of_par=sf_charalloc(strlen(par));
  strcpy (copy_of_par, par);
  list_of_keys[0]=strtok(copy_of_par,",");
  for(ikey=1; ikey<*numkeys; ikey++){
    list_of_keys[ikey]=strtok(NULL,",");
  }
  /*should:
    free(copy_of_par);
    free(par);
    but it causes part of list_of_keys to be clobberred 
    this is a very minor memory leak
  */
  return list_of_keys;
}

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  int n_traces; 
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  int* iheader=NULL;
  int i_trace;
  int tempint;
  char type_input_record[5];
  int input_record_length;
  int numkeys;
  int ikey;
  char** list_of_keys;
  int *indx_of_keys;
  
  sf_init (argc,argv);

   /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
  fprintf(stderr,"verbose=%d\n",verbose);
 
  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/

  fprintf(stderr,"read in file name\n");  
  in = sf_input ("in");

  fprintf(stderr,"read out file name\n");
  out = sf_output ("out");

  if (!sf_histint(in,"n1_traces",&n1_traces))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n1_headers",&n1_headers)) 
    sf_error("input data does not define n1_headers");

  n_traces=sf_leftsize(in,1);

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  if (SF_INT == typehead) iheader = sf_intalloc(n1_headers);
  else                    fheader = sf_floatalloc(n1_headers);
 
  fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  fprintf(stderr,"call list of keys\n");
  list_of_keys=sf_getstringarray("key",&numkeys);
  if(list_of_keys==NULL)
    sf_error("The required parameter \"key\" was not found.");

  /* print the list of keys */
  for(ikey=0; ikey<numkeys; ikey++){
    fprintf(stderr,"list_of_keys[%d]=%s\n",ikey,list_of_keys[ikey]);
  }
  
/* maybe I should add some validation that n1== n1_traces+n1_headers+2
     and the record length read in the second word is consistent with 
     n1.  */

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /* the list header keys (including any extras) */
  segy_init(n1_headers,in);
  indx_of_keys=sf_intalloc(numkeys);
  for (ikey=0; ikey<numkeys; ikey++){
    /* kls need to check each of these key names are in the segy header and
       make error message for invalid keys.  Of does segykey do this? NO, just
       segmentation fault. */
    indx_of_keys[ikey]=segykey(list_of_keys[ikey]);
  }



  fprintf(stderr,"start trace loop\n");
  for (i_trace=0; i_trace<n_traces; i_trace++){

    if(i_trace<5)fprintf(stderr,"i_trace=%d, input_record_length=%d\n",
			         i_trace   , input_record_length);
    /* kls. Factor this into a function sf_get_trace */
    sf_charread(type_input_record,4,in);
    type_input_record[4]='\0';
    if(verbose>1)fprintf(stderr,"type_input_record=%s\n",
	                 type_input_record);
    sf_intread(&input_record_length,1,in);
    if(strcmp(type_input_record,"tah ")!=0){
      /* not my kind of record.  Just write it back out */
      sf_error("non tah record found. Better write this code segment\n");
    } else {
      /* process the tah. */
      /* this program just prints selected header keys */

      /**************************/
      /* read trace and headers */
      /**************************/
      if(verbose>1)fprintf(stderr,"read header\n");
      if(SF_INT == typehead) sf_intread  (iheader,n1_headers,in);
      else                   sf_floatread(fheader,n1_headers,in);
      if(verbose>1)fprintf(stderr,"read intrace\n");
      sf_floatread(intrace,n1_traces,in);
      
      if(verbose>1)fprintf(stderr,"tah read.  process and write it.\n");

      for(ikey=0; ikey<numkeys; ikey++){
	fprintf(stderr," %s=",list_of_keys[ikey]);
	if(typehead == SF_INT)fprintf(stderr,"%d",iheader[indx_of_keys[ikey]]);
	else                  fprintf(stderr,"%f",fheader[indx_of_keys[ikey]]);
      }
      fprintf(stderr,"\n");

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
  }

  exit(0);
}

  
