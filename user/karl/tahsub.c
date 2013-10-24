#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../api/c/file.h"
#include <rsf.h>
#include <rpc/xdr.h>


extern int verbose;


/* why do I need to:
  #include "../../api/c/file.h"
  seems like using <> instead of "" should be able to find code.
*/

/* seems like include file.h should take care of these typedef and struct */
#ifndef _sf_file_h

#define SF_MAX_DIM 9
/*^*/

typedef struct sf_File *sf_file;
/*^*/

typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX, SF_SHORT, SF_DOUBLE} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;
/*^*/

#endif

struct sf_File {
    FILE *stream, *head; 
    char *dataname, *buf, *headname;
    sf_simtab pars;
    XDR xdr;
    enum xdr_op op;
    sf_datatype type;
    sf_dataform form;
    bool pipe, rw, dryrun;
};








/* collection of functions shared by the Mtah programs */

/* global variables. need a temp array for sf_get_tah to pass non tah
   records read from pipe.  This means sf_get_tah is not thread safe,
   but the intended use is to read stdin, so why whould you want to 
   run multiple sf_get_tah on more then one thread? */
char* sf_get_tah_bytebuffer=NULL;
int sf_get_tah_bytebuffer_len=0;


/**********************************************************************/
int sf_getnumpars(const char* key){
  char* par;
  char* one_par;
  char* copy_of_par;
  int numpars=0;
  
  fprintf(stderr,"in sf_getnumpars call sf_getstring key=%s\n",key);
  par=sf_getstring(key);
  fprintf(stderr,"back from sf_getstring par=%s\n",par);
  /* this line may not do anything but make code easier to understand */
  if(par==NULL) return 0; /* if key is not an input parameter return 0 */
   
  copy_of_par=sf_charalloc(strlen(par));
  strcpy (copy_of_par, par); /* copy par, strtok destroys it */
  fprintf(stderr,"copy_of_par=%s\n",copy_of_par);
  
  one_par=strtok(copy_of_par,",");
  while (one_par!=NULL){
    fprintf(stderr,"before incrmeent numpar=%d\n",numpars); 
    numpars++;
    one_par=strtok(NULL,",");
  }
  /* kls this may clobber something */
  fprintf(stderr,"free(copy_of_par)\n");
  free(copy_of_par);
  fprintf(stderr,"free(par)\n");
  free(par);
  
  return numpars;
}
 
  
char** sf_getnstring(char* key, int* numkeys){
  /* perhaps this function should go with functions like sf_getint 
     and sf_getstring.  For now I'll keep here with the tah_io programs */

  /* this function reads parameters from the command line.  Original 
     application was in sftahgethw where you used it to read a list
     of headers keys eg key=ep,fldr,offset,sx,sy,gx,gy */
  /* kls decompose this into two functions.  
     I think I can use sf_getstrings :
     bool sf_getstrings (const char* key,@out@ char** par,size_t n) 
     I think par needs to be allocated before calling.  just like sf_getfloats:
     bool sf_getfloats (const char* key,@out@ float* par,size_t n) 
       n, then number of comma seperated values needs to be computed so
       the array par can be allocated.  maybe a new function sf_getnumkeys 
  */
         
  char** list_of_keys;
  char* par;
  char* one_par;
  char* copy_of_par;
  int ikey; 
  
  par=sf_getstring(key);

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
     list_of_keys. this spooked me on exactly how strtok works */ 
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

/**********************************************************************/
int sf_put_tah(float* trace, float* header, 
	       int n1_traces, int n1_headers, sf_file file){
  int input_record_length;
  
  fwrite("tah ",sizeof(char),4,file->stream);

  input_record_length=sizeof(float)*(n1_traces+n1_headers);
  /*  fprintf(stderr,"sf_put_tah write input_record_length=%d\n",
      input_record_length); */

  fwrite(&input_record_length,sizeof(int),1,file->stream);


  if(n1_traces != fwrite(trace,sizeof(float), n1_traces, file->stream))
    return 1;

  if(n1_headers != fwrite(header,sizeof(float), n1_headers, file->stream))
    return 1;
  
  return 0;
}

/**********************************************************************/
int sf_get_tah(float* trace, float* header, 
	       int n1_traces, int n1_headers, sf_file file){
  char type_input_record[5];
  int input_record_length;

  if(4!=fread(type_input_record,sizeof(char),4,file->stream))
    return 1;
  type_input_record[4]='\0';
  /* fprintf(stderr,"sf_get_tah type_input_record=\"%s\"\n",
     type_input_record); */
  if(1!=fread(&input_record_length,sizeof(int),1,file->stream))
    return 1;
  /*   fprintf(stderr,"sf_put_tah write input_record_length=%d\n",
       input_record_length); */
  if(strcmp(type_input_record,"tah ")!=0){
    /* not my kind of record.  Just write it back out */
    /* this is where I need to add the sf_get_tah_buffer allocations kls */
    sf_error("non tah record found. Better write this code segment\n");
  }
  
  if(sizeof(float)*(n1_traces+n1_headers)!=input_record_length){
    fprintf(stderr,"in sf_get_tah.  n1_traces, and n1_headers are\n");
    fprintf(stderr,"inconsistent with input record length\n");
    fprintf(stderr,"n1_traces=%d, n1_headers=%d, sizeof(float)=%d\n",
	    n1_traces   , n1_headers   , (int)sizeof(float) );
    fprintf(stderr,"input_record_length=%d\n",input_record_length);
    sf_error("sizeof(float)*(n1_traces+n1_headers)!=input_record_length\n");
  }
  
  /*fprintf(stderr,"read n1_traces=%d floats from input file\n",n1_traces); */
  if(n1_traces!=fread(trace,sizeof(float), n1_traces, file->stream))
    return 1;

  /*fprintf(stderr,"read n1_headers=%d floats from input file\n",n1_headers);*/
  if(n1_headers!=fread(header,sizeof(float), n1_headers, file->stream))
    return 1;
  /*fprintf(stderr,"return from sf_get_tah\n"); */
  return 0;
}



void sf_tahwritemapped(float* trace, void* iheader, 
		       int n1_traces, int n1_headers,
		       sf_file output,sf_file outheaders,
		       sf_datatype typehead,sf_axis* output_axa_array,
		       int* indx_of_keys, int dim_output,
		       off_t* n_output, off_t* n_outheaders){

  int iaxis;
  double header[SF_MAX_DIM];
  off_t i_output[SF_MAX_DIM];
  off_t file_offset;
  bool trace_fits_in_output_file;

  if(verbose>2)
    for (iaxis=0; iaxis<SF_MAX_DIM; iaxis++){
    fprintf(stderr,"axis=%d sf_n(output_axa_array[iaxis])=%d\n",
 	    iaxis+1,sf_n(output_axa_array[iaxis]));
  }

  /* compute the i_output, the output trace index array.  Use the trace 
     headers and axis definitions.  Axis definitions are label#, n#, o#, 
     and d# */
  /* make sure the trace fits in the output file.  A trace does not fit if
     if falls outside the range defined by n#, o#, and d#.  A trace must fall
     on the o#, d# locations (ie if you say o1=10 d1=10 n1=10, traces -500, 
     10000, and 15 do not fit in the output file. */

  trace_fits_in_output_file=true;

  i_output[0]=0;
  for (iaxis=1; iaxis<dim_output; iaxis++){
    double doubleindex;
    if(SF_INT == typehead) 
      header[iaxis]=(double)(((int  *)iheader)[indx_of_keys[iaxis]]);
    else                   
      header[iaxis]=(double)(((float*)iheader)[indx_of_keys[iaxis]]);
    doubleindex=(header[iaxis]-sf_o(output_axa_array[iaxis]))/
                sf_d(output_axa_array[iaxis]);
    i_output[iaxis]=llround(doubleindex);
    if(0.01<fabs(doubleindex-i_output[iaxis]) ||
       i_output[iaxis]<0 ||
       i_output[iaxis]>=sf_n(output_axa_array[iaxis]) ){
	trace_fits_in_output_file=false;
	fprintf(stderr,"trace rejected\n");
	break;
    }
    if(verbose>2){
      fprintf(stderr,"iaxis=%d n=%d o=%f d=%f header=%f \n",
	      iaxis,
	      sf_n(output_axa_array[iaxis]),
	      sf_o(output_axa_array[iaxis]),
	      sf_d(output_axa_array[iaxis]),
	      header[iaxis]);
    }
  }
  /* kls need to check doubleindex is about the same as i_output[iaxis]
     this will allow you to write every third trace to an output file
     by using a large increment.  This is good to write gathers every km 
     for velocity analysis of qc.  
     kls also need to check index is in the range 
     [0,sf_n(output_axa_array[iaxis])
  */
    /* kls why does this fail?
      fprintf(stderr,"axis=%d output_axa_array[iaxis]->n=%d\n",
	  iaxis+1,output_axa_array[iaxis]->n); */
  if(trace_fits_in_output_file){
    file_offset=sf_large_cart2line(dim_output,n_output,i_output)*sizeof(float);
    if(verbose>2)fprintf(stderr,"file_offset=%lld\n",(long long) file_offset);
    sf_seek(output,file_offset,SEEK_SET);
    sf_floatwrite(trace,n1_traces,output);
    
    file_offset=sf_large_cart2line(dim_output,n_outheaders,i_output)*
                sizeof(float);
    sf_seek(outheaders,file_offset,SEEK_SET);
    if(SF_INT == typehead) sf_intwrite  ((int  *)iheader,n1_headers,outheaders);
    else                   sf_floatwrite((float*)iheader,n1_headers,outheaders);
    if(verbose>2){
      for(iaxis=2; iaxis<dim_output; iaxis++){
	fprintf(stderr,"indx_of_keys[%d]=%d ",
		iaxis,indx_of_keys[iaxis]);
	if(typehead == SF_INT)fprintf(stderr,"%d\n",
				      ((int  *)iheader)[indx_of_keys[iaxis]]);
	else                  fprintf(stderr,"%g\n",
				      ((float*)iheader)[indx_of_keys[iaxis]]);
      }
      fprintf(stderr,"\n");
    }
  }
  
}
