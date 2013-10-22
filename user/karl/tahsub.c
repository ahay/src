/* collection of functions shared by the Mtah programs */

/* global variables. need a temp array for sf_get_tah to pass non tah
   records read from pipe.  This means sf_get_tah is not thread safe,
   but the intended use is to read stdin, so why whould you want to 
   run multiple sf_get_tah on more then one thread? */
char* sf_get_tah_bytebuffer=NULL;
int sf_get_tah_bytebuffer_len=0;

char** sf_getstringarray(char* key, int* numkeys){
  /* perhaps this function should go with functions like sf_getint 
     and sf_getstring.  For now I'' keep if here with the tah_io programs */

  /* this function reads parameters from the command line.  Original 
     application was in sftahgethw where you used it to read a list
     of headers keys eg key=ep,fldr,offset,sx,sy,gx,gy */
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
int sf_get_tah(float* trace, float* header, 
	       int n1_traces, int n1_headers, sf_file file){
  char type_input_record[5];
  
  fread(type_input_record,sizeof(char),4,file->stream);
  type_input_record[4]='\0';
  if(verbose>1)fprintf(stderr,"type_input_record=%s\n",
		       type_input_record);
  sf_intread(&input_record_length,1,in);
  if(strcmp(type_input_record,"tah ")!=0){
    /* not my kind of record.  Just write it back out */
    /* this is where I need to add the sf_get_tah_buffer allocations kls */
    sf_error("non tah record found. Better write this code segment\n");
  }

zas






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
} /* end of block that is to be moved into sf_get_tah function */


if(


