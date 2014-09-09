#ifndef __GRIDPP_CREATE__
#define __GRIDPP_CREATE__

#include "usempi.h"
#include "iwave_fopen.h"
#include "except.hh"
#include "std_cpp_includes.hh"

using namespace RVL;
using namespace std;

// prec=0: float
// prec!=0: double

void create_hfile(string hfile,int prec=0) {

  if (hfile.size()>120) {
    RVLException e;
    e<<"Error: create_hfile\n";
    e<<"filename "<<hfile<<" longer than 120 chars\n";
    throw e;
  }

  char * fname=(char *)malloc(128*sizeof(char));
  FILE * fp;
  string dfile = hfile+"@";

  // first try to open for read
  strcpy(fname,hfile.c_str());
  if (fp = iwave_fopen(&fname,"r",NULL,stderr)) {
    strcpy(fname,dfile.c_str());
    fp = iwave_fopen(&fname,"r",NULL,stderr);
  }

  // if either header or data file not present, create both

  if (!fp) {

    strcpy(fname,hfile.c_str());
    
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file testgrid.rsf not opened\n";
      throw e;
    }

    int n1=11;
    int n2=21;
    int n3=41;
    float d1=0.1;
    float d2=0.1;
    float d3=1.0;
    float o1=-0.5;
    float o2=-1.0;
    float o3=0.0;

    int dim = 2;
    int gdim = 3;

    fprintf(fp,"n1=%d\n",n1);
    fprintf(fp,"d1=%e\n",d1);
    fprintf(fp,"o1=%e\n",o1);
      
    fprintf(fp,"n2=%d\n",n2);
    fprintf(fp,"d2=%e\n",d2);
    fprintf(fp,"o2=%e\n",o2);
      
    fprintf(fp,"n3=%d\n",n3);
    fprintf(fp,"d3=%e\n",d3);
    fprintf(fp,"o3=%e\n",o3);

    fprintf(fp,"dim=%d\n",dim);
    fprintf(fp,"gdim=%d\n",gdim);

    fprintf(fp,"data_format=native_float\n");

    if (prec) fprintf(fp,"data_format=native_double\n");
    else fprintf(fp,"data_format=native_float\n");

    fprintf(fp,"data_type = fungus\n");
    fprintf(fp,"in=%s\n",dfile.c_str());
      
    iwave_fclose(fp);
    fflush(fp);

    strcpy(fname,dfile.c_str());
    float * buf = NULL;
    double * dbuf = NULL;
    if (prec) {
      dbuf=(double *)malloc(n1*n2*n3*sizeof(double));
      for (int i=0;i<n1*n2*n3;i++) dbuf[i]=0.0;
    }
    else {
      buf=(float *)malloc(n1*n2*n3*sizeof(float));
      for (int i=0;i<n1*n2*n3;i++) buf[i]=0.0;
    }
      
    //      FILE * 
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file testgrid.rsf not opened\n";
      throw e;
    }
    if (prec) fwrite(dbuf,sizeof(double),n1*n2*n3,fp);
    else fwrite(buf,sizeof(float),n1*n2*n3,fp);

    iwave_fclose(fp);

    fflush(fp);

    // reopen prototype files for read
    strcpy(fname,hfile.c_str());
    fp = iwave_fopen(&fname,"r",NULL,stderr);
    strcpy(fname,dfile.c_str());
    fp = iwave_fopen(&fname,"r",NULL,stderr);

  }

  free(fname);
}

float dot() {
  int n1=11;
  int n2=21;
  int n3=41;
  float d1=0.1;
  float d2=0.1;
  float d3=0.1;
  return n1*n2*n3*d1*d2*d3;
}

#endif
