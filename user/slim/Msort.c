/* Sort a float/complex vector by absolute values.

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

static int asc(const void * a, const void * b)
{
  return (*(float*)a<*(float*)b)? -1: (*(float*)a>*(float*)b)? 1:0;
}

static int desc(const void * a, const void * b)
{
  return (*(float*)a>*(float*)b)? -1: (*(float*)a<*(float*)b)? 1:0;
}

static float *dataloader(sf_file in, int n, bool complex_data)
{
    int i;
    float *data;
    sf_complex *cdata;
    
    data = sf_floatalloc(n); 
    if (complex_data){
	/* read in memory and compute absolute value of vector entries */
	cdata = sf_complexalloc(n);
	sf_complexread(cdata,n,in);
	for (i=0; i<n; i++){
	    data[i] = cabsf(cdata[i]);
	}
	free(cdata);
    } else{
	/* read in memory and compute absolute value of vector entres */
	sf_floatread(data,n,in);
	for (i=0; i<n; i++){
	    data[i] = fabsf(data[i]);
	}
    }
  
    return(data);
}

int main(int argc, char* argv[])
{
  int n, nbchunks, nfit, i, j;
  off_t n1, memsize;
  int mem; /* for avoiding int to off_t typecast warning */
  float *data;
  float *tmp;
  FILE **fp;
  char **fn;
  int nloop;
  float *currentext;
  int icurrentext;
  bool complex_data = false; /* true if input data is complex */
  bool ascmode;
  sf_file in, out; /* Input and output files */

  /* Initialize RSF */
  sf_init(argc,argv);

  /* standard input */
  in = sf_input("in");
  /* standard output */
  out = sf_output("out");

  if (!sf_getint("memsize",&mem))
      mem=sf_memsize();
  /* Max amount of RAM (in Mb) to be used */
  memsize = mem * (1<<20); /* convert Mb to bytes */

  if (!sf_getbool("ascmode",&ascmode)) ascmode = false;
  /* y=ascending; n=descending */

  /* n1 is the vector length */
  n1 = sf_filesize(in);
  
  /* check that the input is float or complex and format output appropriately */
  if (SF_COMPLEX == sf_gettype(in)){
      complex_data = true;
      sf_settype(out,SF_FLOAT);
  } else if (SF_FLOAT != sf_gettype(in)) 
      sf_error("Need float or complex input");
  
  if (n1*sizeof(float)>memsize){
    /* OUT-OF-CORE MODE */
    sf_warning("Going out of core...\n"
	       "Increase memsize for in-core");
    
    n = memsize/sizeof(float); /* # of element per chunk */
    nbchunks = n1/n+1; /* # of chunks */
    nfit = n1%n; /* # of element in the last chunk */
    
    fp = (FILE**) sf_alloc(nbchunks,sizeof(FILE *));
    fn = (char**) sf_alloc(nbchunks,sizeof(char *));

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
      
      fp[i] = sf_tempfile(&fn[i],"w+b");
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
    
    currentext = sf_floatalloc(1);
    for (i=0; i<n1; i++){
      *currentext = data[0];
      icurrentext = 0;
      /* find current extremum */
      for (j=0; j<nbchunks; j++){
	if (ascmode){
	  if (data[j]<*currentext){
	    *currentext = data[j];
	    icurrentext = j;
	  }
	}
	else{
	  if (data[j]>*currentext){
	    *currentext = data[j];
	    icurrentext = j;
	  }
	}
      }
      
      fread(tmp, sizeof(float),1, fp[icurrentext]);
      if (!(feof(fp[icurrentext]))){
	data[icurrentext] = *tmp;
      }
      else{
	if (ascmode)
	  data[icurrentext] = FLT_MAX;
	else
	  data[icurrentext] = 0;
      }
      
      sf_floatwrite(currentext,1,out);
    }
    
    for (i=0; i<nbchunks; i++) {
      fclose(fp[i]);
      unlink(fn[i]);
    }
  }
  else{
    /* IN-CORE MODE */
    data = dataloader(in, n1, complex_data);
    
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
  
