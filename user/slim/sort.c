/*
Sort in ascending/descending order absolute value entries of a
float/complex vector.

Written by: Gilles Hennenfent & Henryk Modzelewski, UBC
Created: February 2006
*/

/*
Copyright (C) 2006 The University of British Columbia at Vancouver

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

int asc(const void * a, const void * b)
{
  return (*(float*)a<*(float*)b)? -1: (*(float*)a>*(float*)b)? 1:0;
}

int desc(const void * a, const void * b)
{
  return (*(float*)a>*(float*)b)? -1: (*(float*)a<*(float*)b)? 1:0;
}

float *dataloader(sf_file in, int n, bool complex_data)
{
  int i;
  float *data;
  data = sf_floatalloc(n); 
  if (complex_data){
    /* read in memory and compute absolute value of vector entries */
    float complex *cdata;
    cdata = sf_complexalloc(1);
    for (i=0; i<n; i++){
      sf_complexread(cdata,1,in);
      data[i] = sqrt( pow(creal(cdata[0]),2) + pow(cimag(cdata[0]),2));
    }
  }
  else{
    /* read in memory and compute absolute value of vector entres */
    sf_floatread(data,n,in);
    for (i=0; i<n; i++){
      data[i] = sqrt(pow(data[i],2));
      }
  }
  
  return(data);
}

int main(int argc, char* argv[])
{
  int n1, mem, memsize, n, nbchunks, nfit, i, j;
  size_t len;
  float *data;
  float *tmp;
  FILE **fp;
  int nloop;
  char *rformat, *cformat;
  float *currentmax;
  int icurrentmax = 0;
  bool complex_data = false;// true if input data is complex
  bool ascmode;// true if input data is complex
  sf_file in, out; /* Input and output files */
  
  /* Initialize RSF */
  sf_init(argc,argv);
  
  /* standard input */
  in = sf_input("in");
  /* standard output */
  out = sf_output("out");
  
  if (!sf_getint("memsize",&mem)) mem = 500;
  /* Available memory size (in Mb) */
  memsize = mem * (1 << 20); /* convert Mb to bytes */
  
  if (!sf_getbool("ascmode",&ascmode)) ascmode = false;
  /* y=ascending; n=descending */

  /* n1 is the vector length */
  if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
  /* leftsize gets n2*n3*n4*... */
  if (sf_leftsize(in,1) != 1) sf_error("Need a vector as input");    
  
  /* check that the input is float or complex and format output appropriately */
  if (SF_COMPLEX == sf_gettype(in)){
    complex_data = true;
    cformat = sf_histstring(in,"data_format");
    len = strlen(cformat)+1;
    rformat = sf_charalloc(len);
    memcpy(rformat,cformat,len);
    strcpy(strstr(rformat,"complex"),"float");
    sf_setformat(out,rformat);
  }
  else if (SF_FLOAT != sf_gettype(in)) sf_error("Need float or complex input");
  
  if (n1*sizeof(float)>memsize){
    /* OUT-OF-CORE MODE */
    sf_warning("Going out of core...\n"
	       "Increase memsize for in-core");
    
    n = memsize/sizeof(float); /* # of element per chunk */
    nbchunks = n1/n+1; /* # of chunks */
    nfit = n1%n; /* # of element in the last chunk */
    
    /* read chunk, sort and write in temp files */
    if (!getenv("TMPDIR"))
      {
	putenv("TMPDIR=./");
	sf_warning("Set TMPDIR to ./");
      }
    else{
      sf_warning("Use user-defined TMPDIR = %s",getenv("TMPDIR"));
    }


    fp = malloc(nbchunks*sizeof(FILE *));
    
    for (i=0; i<nbchunks; i++){
      if (i == nbchunks-1) nloop = nfit;
      else nloop= n;
      
      data = dataloader(in, nloop, complex_data);
      
      /* sort data */
      if (ascmode){
	qsort(data,nloop,sizeof(float),asc);
      }
      else{
	qsort(data,nloop,sizeof(float),desc);
      }
      
      fp[i] = tmpfile();
      fwrite(data, sizeof(float), nloop, fp[i]);
      rewind(fp[i]);
    }
    
    /* merge sorted vectors */
    tmp = sf_floatalloc(1);
    data = sf_floatalloc(nbchunks);
    for (i=0; i<nbchunks; i++){
      fread(tmp, sizeof(float), 1, fp[i]);
      data[i] = *tmp;
    }
    
    currentmax = sf_floatalloc(1);
    for (i=0; i<n1; i++){
      *currentmax = 0;
      /* find current max */
      for (j=0; j<nbchunks; j++){
	if (data[j]>*currentmax){
	  *currentmax = data[j];
	  icurrentmax = j;
	}
      }
      
      fread(tmp, sizeof(float),1, fp[icurrentmax]);
      if (!(feof(fp[icurrentmax]))){
	data[icurrentmax] = *tmp;
      }
      else{
	data[icurrentmax] = 0;
      }
      
      sf_floatwrite(currentmax,1,out);
    }
    
    for (i=0; i<nbchunks; i++)
      fclose(fp[i]);
    
  }
  else{
    /* IN-CORE MODE */
    float *data;
    data = dataloader(in, n1, complex_data);
    sf_fileclose(in);
    
    /* sort data */
    if (ascmode){
      qsort(data,n1,sizeof(float),asc);
    }
    else{
      qsort(data,n1,sizeof(float),desc);
    }
    
    /* write sorted data to file */
    sf_floatwrite(data,n1,out);
  }  

  sf_close();
  
  exit(0);
}
  
