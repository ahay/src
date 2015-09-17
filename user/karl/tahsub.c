/* collection of functions shared by the Mtah programs */
/*
  Copyright (C) 2013 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include "tahsub.h"

/* global variables. need a temp array for sf_get_tah to pass non tah
   records read from pipe.  This means sf_get_tah is not thread safe,
   but the intended use is to read stdin, so why whould you want to 
   run multiple sf_get_tah on more then one thread? 
char* sf_get_tah_bytebuffer=NULL;
int sf_get_tah_bytebuffer_len=0;
*/
int counttokens(const char* delimittedtokens, const char* delimiters)
/*< return number of delimiter seperated tokens in string>*/
{
  char* one_token;
  char* copy_of_delimittedtokens;
  int numtoken=0;

  if(NULL==delimittedtokens)return 0;
  /* allocate and copy delimittedtokens
     strtok destroys the string as it splits it up */
  copy_of_delimittedtokens=sf_charalloc(strlen(delimittedtokens)+1);
  strcpy (copy_of_delimittedtokens, delimittedtokens);
  fprintf(stderr,"copy_of_delimittedtokens=%s\n",copy_of_delimittedtokens);
  
  one_token=strtok(copy_of_delimittedtokens,delimiters);
  while (one_token!=NULL){
    if(false)fprintf(stderr,"before increment numtoken=%d\n",numtoken); 
    numtoken++;
    one_token=strtok(NULL,delimiters);
  }
  if(false)fprintf(stderr,"free(copy_of_delimittedtokens)\n");
  /* did not save and of the pointers returned by strtok, so OK to free! */
  free(copy_of_delimittedtokens);
  if(false)fprintf(stderr,"return from counttokens numtoken=%d\n",numtoken);
  return numtoken;
}

char** delimittedtokens2list(const char* delimittedtokens, 
			     const char* delimiters, 
			     int* numtoken)
/*<return number of delimitter seperated tokens in string and list of tokens>*/
{
  int itoken;
  char** list_of_tokens;
  char* copy_of_delimittedtokens;
  
  /* count number of keys */
  *numtoken=counttokens(delimittedtokens,delimiters);
  if(false)fprintf(stderr,"in delimittedtokens2list *numtoken=%d\n",*numtoken);
  if(0==*numtoken)return NULL;
  /* break the long string into an array  of strings */
  list_of_tokens=(char**)sf_alloc(*numtoken,sizeof(char*)); 
  /* allocate and copy delimittedtokens
     strtok destroys the string as it splits it up */
  copy_of_delimittedtokens=sf_charalloc(strlen(delimittedtokens)+1);
  strcpy (copy_of_delimittedtokens, delimittedtokens);

  list_of_tokens[0]=strtok(copy_of_delimittedtokens,delimiters);
  for(itoken=1; itoken<*numtoken; itoken++){
    list_of_tokens[itoken]=strtok(NULL,delimiters);
  }
  /* 
     I think list_of_tokens points into copy_of_delimittedtokens. 
     Do not free!
  */
  return list_of_tokens;  
}

int getnumpars(const char* key)
/*< get number of pars >*/
{
  char* par;
  int numpars=0;
  
  sf_warning("in getnumpars call sf_getstring key=%s\n",key);
  if(NULL==(par=sf_getstring(key)))return 0;
   
  numpars=counttokens(par,",");

  fprintf(stderr,"free(par)\n");
  free(par);
  
  return numpars;
}

char** sf_getnstring(char* key, int* numkeys)
/*< get strings >*/
{
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
         
  char* par;
    
  if(!(par=sf_getstring(key)))return NULL;

  return delimittedtokens2list(par,",",numkeys);
}

void put_tah(float* trace, float* header, 
	     int n1_traces, int n1_headers, sf_file file)
/*< put tah >*/
{
  int input_record_length;
  
  sf_charwrite("tah ",4,file);

  input_record_length=sizeof(int)*(n1_traces+n1_headers);
  /*  fprintf(stderr,"sf_put_tah write input_record_length=%d\n",
      input_record_length); */

  sf_intwrite(&input_record_length,1,file);
  sf_floatwrite(trace,n1_traces,file);
  sf_floatwrite(header,n1_headers,file);
}

bool get_tah(float* trace, float* header, 
	     int n1_traces, int n1_headers, sf_file file)
/*< get tah.  return 1 if eof encountered, 0 otherwise >*/
{
  int input_record_length;
  char type_input_record[5];

  /* sergey suggests
     if (!sf_try_charread("tah",file)) return 1;
     but I cannot find function and I donot think it will do all that I need
     For now I am changing back to my stub kls */
  
  if(4!=sf_try_charread2(type_input_record,4,file)){
    /* must have encounterred eof. */
    return true;
  }
  type_input_record[4]='\0';
  sf_intread(&input_record_length,1,file);
  /* fprintf(stderr,"sf_get_tah type_input_record=\"%s\"\n",
     type_input_record); */
  if(strcmp(type_input_record,"tah ")!=0){
    /* not my kind of record.  Just write it back out */
    /* Right now tah programs only have tah records.  Other type records 
       used to be be helpful to me.  They allowed me to share the pipe.  
       Program looked for the type records and they wanted and ignored the 
       ones it was not interested in an may did not even know.  If I 
       add other record types this is where I need to add the 
       sf_get_tah_buffer allocations kls */
    sf_error("non tah record found. Better write this code segment\n");
  }     
  
  if(sizeof(float)*(n1_traces+n1_headers)!=input_record_length){
      sf_warning("%s: n1_traces, and n1_headers are",__FILE__);
      sf_warning("inconsistent with input record length");
      sf_warning("n1_traces=%d, n1_headers=%d, sizeof(float)=%d",
		  n1_traces   , n1_headers  ,(int)sizeof(float) );
      sf_warning("input_record_length=%d",input_record_length);
      sf_error("sizeof(float)*(n1_traces+n1_headers)!=input_record_length");
  }
  
  /*fprintf(stderr,"read n1_traces=%d floats from input file\n",n1_traces); */

  sf_floatread(trace,n1_traces,file);

  /*fprintf(stderr,"read n1_headers=%d floats from input file\n",n1_headers);*/

  sf_floatread(header,n1_headers,file);

  return false;
}


void tahwritemapped(int verbose, float* trace, void* iheader, 
		    int n1_traces, int n1_headers,
		    sf_file output,sf_file outheaders,
		    sf_datatype typehead,
		    sf_axis* output_axa_array,
		    int* indx_of_keys, int dim_output,
		    off_t* n_output, off_t* n_outheaders)
/*< tah write mapped >*/
{

  int iaxis;
  double header[SF_MAX_DIM];
  off_t i_output[SF_MAX_DIM];
  off_t file_offset;
  bool trace_fits_in_output_file;
  static int number_traces_rejected=0;

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
	number_traces_rejected++;
	if(number_traces_rejected<10){
	  fprintf(stderr,"trace rejection #%d\n",number_traces_rejected);
	}
	if(number_traces_rejected==10){
	  fprintf(stderr,"another rejected trace.  reporting will stop.\n");
	}
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

void tahwriteseq(int verbose, float* trace, void* iheader, 
		 int n1_traces, int n1_headers,
		 sf_file output,sf_file outheaders,
		 sf_datatype typehead)
/*< tah write seq >*/
{

  /* write trace and header to the output files.  LIke tahwritemapped, but  
     no seeks;  just write to the next location */

  sf_floatwrite(trace,n1_traces,output);
  
  if(SF_INT == typehead) sf_intwrite  ((int  *)iheader,n1_headers,outheaders);
  else                   sf_floatwrite((float*)iheader,n1_headers,outheaders);
  if(verbose>2){
      fprintf(stderr,"trace and header written to output file\n");
  }
}

