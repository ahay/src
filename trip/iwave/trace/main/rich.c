#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

#define TOL 0.001

char * sdoc[] = { 
  "Usage: richext.x in1= in2= out= par1= par2= order=[1]",
  "",
  "Synopsis: assumes that SU trace data files in1 and in2 represent",
  "same parameter-dependent quantities, with two different values of the",
  "parameter given by par1 and par2. Assumes that these inputs asymptotically",
  "approximate the same limit quantities, with convergence rate given",
  "by the input order.Uses Richardson extrapolation to estimate the error",
  "in the 2nd input in2, and outputs extrapolated estimate of limit to file",
  "out. Also prints out a list of estimated mean-square errors for in2 for",
  "each trace, and a mean-square error estimate the entire in2 file. Since",
  "errors are computed for in2, typically par1 > par2 so that in2 is the",
  "presumptively more accurate input."
  "",
  "The numbers of time samples in each corresponding trace must be the same,",
  "and the time sample rates are checked to be the same to within 0.1%. No",
  "other headers are checked for compatibility between the two data sets.",
  "Traces are compared until one file or the other runs out of traces.",
  "",
  "Required inputs:",
  "  in1 [string]  = data computed for parameter value par1",
  "  in2 [string]  = data computed for parameter value par2",
  "  out [string]  = Richardson-extrapolated output traces",
  "  par1 [float]  = parameter for in1",
  "  par2 [float]  = parameter for in2",
  "  order  [int]  = assumed order of convergence (default = 1)",
  NULL};

int main(int argc, char ** argv) {
  
  PARARRAY * par;    /* param array */
 
  char * in1;      /* input 1 file name */
  char * in2;      /* input 2 file name */
  char * out;      /* output file name */
  FILE * fp1;      /* input 1 file pointer */
  FILE * fp2;      /* input 2 file pointer */
  FILE * fp;       /* output file pointer */
  segy tr1;        /* input 1 trace workspace */
  segy tr2;        /* input 2 trace workspace */
  segy tr;         /* output trace workspace */
  float h1;        /* input 1 param */
  float h2;        /* input 2 param */
  float h2p;       /* (input 2 param)^order */
  Value val;       /* header word workspace */
  float dt=0.0f;   /* time step */
  int nt;          /* number of samples per trace */
  int order=1;     /* asymptotic order of data */
  double wk;       /* double workspace */
  float rfac;      /* Richardson factor: 1/(h1^order-h2^order) */
  float esamp;     /* sample error */
  float gpow=0.0f; /* power of input 2 data */
  float tpow;      /* power of input 2 trace */
  float terror;    /* trace error estimate */
  float gerror;    /* data error estimate */
  int i,itr;       /* counters */

  xargc=argc; xargv=argv;
  requestdoc(1);

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    printf("Error parsing input data. ABORT.\n");
    return 1;
  }
             
  if (ps_flcstring(*par,"in1",&in1)) {
    printf("Error reading in1. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"in2",&in2)) {
    printf("Error reading in2. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"out",&out)) {
    printf("Error reading out. ABORT.\n");
    exit(1);
  }

  if (ps_flfloat(*par,"par1",&h1)) {
    printf("Error reading par1. ABORT.\n");
    exit(1);
  }

  if (ps_flfloat(*par,"par2",&h2)) {
    printf("Error reading par1. ABORT.\n");
    exit(1);
  }

  ps_flint(*par,"order",&order);

  /* compute Richardson factor */
  if (fabs(h1-h2) < iwave_min(TOL*h1,TOL*h2)) {
    printf("Error: par1=%e too close to par2=%e\n",h1,h2);
    exit(1);
  }
  
  wk=h1;
  rfac=pow(wk,order);
  wk=h2;
  h2p=pow(wk,order);
  rfac-=h2p;
  if (rfac < FLT_EPSILON) {
    printf("Error: computed Richardson factor %e too small\n",rfac);
    exit(1);
  }

  rfac=1.0/rfac;

  /* open data files */
  if (!(fp1=fopen(in1,"r"))) {
    printf("Error failed to open 1st input file = %s. ABORT.\n",in1);
    exit(1);
  }

  if (!(fp2=fopen(in2,"r"))) {
    printf("Error failed to open 2nd input file = %s. ABORT.\n",in2);
    exit(1);
  }

  if (!(fp=fopen(out,"w"))) {
    printf("Error failed to open output file = %s. ABORT.\n",out);
    exit(1);
  }

  /* read loop */
  gerror=0.0;
  itr=0;
  printf("******************************************\n");
  printf("* Richardson Mean Square Error Estimator \n");
  printf("*                                        \n");
  printf("* Input 1 = %s, param = %12.6e\n",in1,h1);
  printf("* Input 2 = %s, param = %12.6e\n",in2,h2);
  printf("* Output  = %s\n",out);
  printf("* Assumed convergence order = %d\n",order);
  printf("\n");

  while (fgettr(fp1,&tr1) && fgettr(fp2,&tr2)) {

    /* check that time sampling is compatible */
    gethdval(&tr1,"ns",&val);
    nt = vtoi(hdtype("ns"),val);
    gethdval(&tr2,"ns",&val);
    if (nt != vtoi(hdtype("ns"),val)) {
      printf("Error: mismatch in number of time samples\n");
      exit(1);
    }
    
    gethdval(&tr1,"dt",&val);
    dt = 0.001*vtof(hdtype("dt"),val);
    gethdval(&tr2,"dt",&val);
    if (fabs(dt - 0.001*vtof(hdtype("dt"),val))> TOL*dt) {
      printf("Error: mismatch in time sample rate\n");
      exit(1);
    }

    /* whack in2 header into output trace */
    memcpy(&tr,&tr2,HDRBYTES*sizeof(char));

    /* compute Richardson estimate of mean square difference */
    terror=0.0;
    tpow=0.0;
    for (i=0;i<nt;i++) {
      esamp=rfac*(tr1.data[i]-tr2.data[i])*h2p;
      tr.data[i]=tr2.data[i]-esamp;
      terror+=esamp*esamp;
      tpow+=tr2.data[i]*tr2.data[i];
    } 
    fputtr(fp,&tr);
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
    gpow=gerror/gpow;
  }
  else {
    gpow=FLT_MAX;
  }
  printf("\n");
  printf("est error: abs=%12.6e rel=%12.6e\n",gerror,gpow);
  printf("******************************************\n");

  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  ps_delete(&par);
}
  
