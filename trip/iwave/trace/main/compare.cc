#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

#define DT_TOL 0.001

char * sdoc[] = { 
  "Usage: comp.x in1= in2=",
  "",
  "Purpose: compute and display relative RMS error between two SEGY",
  "data sets, also indicate whether these are the same to specified",
  "precision relative to total energy.",
  "",
  "Required parameters:",
  "  in1 [string]  = filename for first data set",
  "  in2 [string]  = filename for second data set",
  "",
  "Optional Parameters:",
  "  tol [float]   = relative error threshhold for sample-by-sample",
  "                  comparison.",
  "",
  "Synopsis: computes and displays trace-by-trace and whole-file mean",
  "square error between the two input SU trace data files in1 and in2.",
  "Also returns 0 if relative error is less than specified tolerance.",
  "",
  "This utility functions as a rational replacement for sucmp in constructing",
  "regression tests. sucmp returns false (exit(1)) if any pair of ",
  "corresonding samples have a relative difference exceeding threshhold.",
  "Thus floating-point noise changes, such as caused by moving a test suite",
  "to a different platform, can cause a false negative. comp uses for its",
  "base amplitude the maximum absolute amplitude, averaged sample-by-sample",
  "between the two traces. Thus two data sets are regarded as different if",
  "any pair of corresponding samples differ by more than the specified ",
  "tolerance times the maximum sample. This strategy avoids declaring data",
  "sets different due to the presence of trace segments filled with different",
  "floating point noises which invariably have large relative differences.",
  "",
  "In comparing headers, comp behaves exactly as sucmp does. First it checks",
  "that the numbers of time samples in each corresponding trace are the same,",
  "also the time sample rates to within 0.1%. Then the headers are compared",
  "bitwise."
  ""
  "Absolute and relative errors for each trace, and for the data set as a",
  "whole, are printed to stdout. Traces are compared until one file or the",
  "other runs out of traces. On exit, the program returns 0 if the comparison",
  "test described above succeeds, otherwise false."
  "",
  NULL};


int main(int argc, char ** argv) {
  
  /* ARGUMENTS */
  char * in1;      /* input 1 file name */
  char * in2;      /* input 2 file name */
  float tol;       /* tolerance for comparison - default = 10^-6 */

  /* INTERNAL VARIABLES */
  PARARRAY * par;    /* param array */
  FILE * fp1;      /* input 1 file pointer */
  FILE * fp2;      /* input 2 file pointer */
  segy tr1;        /* input 1 trace workspace */
  segy tr2;        /* input 2 trace workspace */
  Value val;       /* header word workspace */
  float dt=0.0f;   /* time step */
  int nt=0;        /* number of samples per trace */
  float esamp;     /* sample error */
  float msamp;     /* max abs */
  float merror;    /* max error */
  float gpow=0.0f; /* global power of input 2 data */
  float tpow;      /* power of input 2 trace */
  float terror;    /* trace error */
  float gerror;    /* global abs error */
  float rerror;    /* global rel error */
  int itr;         /* trace counter */
  int i;           /* sample counter */
  int err;         /* global error flag */

  xargc=argc; xargv=argv;
  requestdoc(1);

  err=0;

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    printf("Error parsing input data. ABORT.\n");
    exit(1);
  }
             
  if (ps_flcstring(*par,"in1",&in1)) {
    printf("COMP reading in1. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"in2",&in2)) {
    printf("COMP reading in2. ABORT.\n");
    exit(1);
  }

  if (ps_flreal(*par,"tol",&tol)) tol=1.e-6;

  /* open data files */
  if (!(fp1=fopen(in1,"r"))) {
    printf("COMP: failed to open 1st input file = %s. ABORT.\n",in1);
    exit(1);
  }

  if (!(fp2=fopen(in2,"r"))) {
    printf("COMP: failed to open 2nd input file = %s. ABORT.\n",in2);
    exit(1);
  }

  /* read loop */
  gerror=0.0;
  itr=0;
  printf("******************************************\n");
  printf("*      COMParison of SEGY data sets      \n");
  printf("*                                        \n");
  printf("* Input 1 = %s\n",in1);
  printf("* Input 2 = %s\n",in2);
  printf("\n");

  msamp=0.0;
  merror=0.0;

  while (fgettr(fp1,&tr1) && fgettr(fp2,&tr2)) {

    /* check that time sampling is compatible */
    gethdval(&tr1,"ns",&val);
    nt = vtoi(hdtype("ns"),val);
    gethdval(&tr2,"ns",&val);
    if (nt != vtoi(hdtype("ns"),val)) {
      printf("COMP: mismatch in number of time samples, trace %d\n",itr);
      err=1;
    }
    
    gethdval(&tr1,"dt",&val);
    dt = 0.001*vtof(hdtype("dt"),val);
    gethdval(&tr2,"dt",&val);
    if (fabs(dt - 0.001*vtof(hdtype("dt"),val))> DT_TOL*dt) {
      printf("COMP: mismatch in time sample rate, trace %d\n",itr);
      err=1;
    }

    if (memcmp(&tr1,&tr2,HDRBYTES)) {
      printf("COMP: headers differ at bit level, trace %d\n",itr);
      err=1;
    }

    if (err) exit(1);

    /* compute mean square difference */
    terror=0.0;
    tpow=0.0;
    for (i=0;i<nt;i++) {
      esamp=tr1.data[i]-tr2.data[i];
      msamp=iwave_max(msamp,0.5*(fabs(tr1.data[i])+fabs(tr2.data[i])));
      merror=iwave_max(merror,fabs(esamp));
      terror+=esamp*esamp;
      tpow+=tr2.data[i]*tr2.data[i];
    } 
    gerror+=terror;
    gpow+=tpow;
    terror=sqrt(dt*terror);
    tpow=sqrt(dt*tpow);
    if (tpow>FLT_EPSILON) {
      tpow=terror/tpow;
    }
    else {
      tpow=FLT_MAX;
    }
    printf("trace %6.3d est error: abs=%12.6e rel=%12.6e\n",itr,terror,tpow);
    itr++;
  }
  gerror=sqrt(dt*gerror);
  gpow=sqrt(dt*gpow);
  if (gpow>FLT_EPSILON) {
    rerror=gerror/gpow;
  }
  else {
    rerror=FLT_MAX;
  }
  printf("\n");
  printf("est error: abs=%12.6e rel=%12.6e\n",gerror,rerror);
  printf("******************************************\n");

  fclose(fp1);
  fclose(fp2);
  ps_delete(&par);

  if (!itr || !nt) exit(1);
  if (merror > tol*msamp) exit(1);
  exit(0);
}
  
