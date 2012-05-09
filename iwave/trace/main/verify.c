/* verify */

#include <verifytrace.h>
#include <par.h>
#include <parser_su.h>

char * sdoc[] = {
  "usage: verify.x trial= comp= [optional parameters]",
  "  trial = filename for test traces",
  "  comp  = filename for comparison traces",
  " ",
  "  optional paramters:",
  "  dim   = [3] dimension - leval values are 2 and 3",
  "  tol   = [0.05] tolerance for normalized max difference",
  "  rc    = [0.1] nominal reflection coefficient, scales direct wave amp",
  "  cref  = [1.5 m/ms] nominal wave velocity",
  "  amp   = [1.0 GPa] nominal amplitude of direct wave at reference distance",
  "  rdist = [1000.0 m] nominal reference distance",
  NULL };

int main(int argc, char ** argv) {

  char * trialfile=NULL;
  char * compfile =NULL;
  int dim=3;
  float tol=0.05;
  float rc=0.1;
  float cref=1.5;
  float amp=1.0;
  float rdist=1000.0;
  float e=0.0;
  float emax=0.0;
  int pf;
  segy trial;
  segy comp;
  FILE * fptrial;
  FILE * fpcomp;
  int ntr;

  ps_initargs(argc,argv);
  xargc=argc; xargv=argv;
  requestdoc(1);

  if (!(ps_getparstring("trial",&trialfile)) ||
      !(ps_getparstring("comp",&compfile))) {
    fprintf(stderr,"Error: testverify\n");
    fprintf(stderr,"missing necessary inputs\n");
    exit(1);
  }

  ps_getparint("dim",&dim);
  ps_getparfloat("tol",&tol);
  ps_getparfloat("rc",&rc);
  ps_getparfloat("cref",&cref);
  ps_getparfloat("amp",&amp);
  ps_getparfloat("rdist",&rdist);

  if (!(fptrial=fopen(trialfile,"r+")) ||
      !(fpcomp =fopen(compfile, "r+"))) {
    fprintf(stderr,"Error: testverify\n");
    fprintf(stderr,"failed to open files\n");
    exit(1);    
  }

  ntr=0;
  emax=0.0;
  while (fgettr(fptrial,&trial) && fgettr(fpcomp,&comp)) {
    pf=verifytrace(&trial,&comp,dim,tol,rc,cref,amp,rdist,&e);
    printf("trace %d: p/f = %d, scaled error = %e\n", 
	   ntr,pf,e);
    emax=fmax(e,emax);
    ntr++;
  }
  printf("max scaled error over all traces = %e\n",emax);
  fclose(fptrial);
  fclose(fpcomp);
}
