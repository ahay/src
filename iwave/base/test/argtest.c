#include <parser.h>

int main(int argc, char ** argv) {

  FILE * fp = NULL;
  PARARRAY * par = ps_new();

  // create fake args
  int xargc = 4;
  char ** xargv = (char **)malloc(xargc*sizeof(char *));
  xargv[0]=NULL;
  xargv[1]=(char *)malloc(strlen("DispFlag=2")+1);
  strcpy(xargv[1],"DispFlag=0");
  xargv[2]=(char *)malloc(strlen("ResTol=0.001")+1);
  strcpy(xargv[2],"ResTol=0.001");
  xargv[3]=(char *)malloc(strlen("par=test/argtest/umin.par")+1);
  strcpy(xargv[3],"par=test/argtest/umin.par");

  // create par file
  char * tname = (char *)malloc(128*sizeof(char));
  memset(tname,'\0',128);
  strcpy(tname,"test/argtest/umin.par");

  if (!(fp = iwave_fopen(&tname,"w",NULL,stderr))) {
    fprintf(stderr,"PANIC!\n");
    exit(1);
  }

  fprintf(fp,"DispFlag=2\n");
  fprintf(fp,"AbsGradTol=0.0\n");
  fprintf(fp,"RelGradTol=0.01\n");
  fprintf(fp,"MaxItn=5\n");
  fprintf(fp,"LS_MaxSample=5\n");
  fprintf(fp,"LS_FirstStep=1.0e-7\n");
  fprintf(fp,"LS_MinStepTol=1.e-10\n");
  fprintf(fp,"LS_FractionOfMaxStep=0.9\n");
  fprintf(fp,"MinDecrease=0.01\n");
  fprintf(fp,"GoodDecrease=0.8\n");
  fprintf(fp,"StepDecrFactor=0.5\n");
  fprintf(fp,"StepIncrFactor=1.8\n");
  fprintf(fp,"BFGS_InvHessianScale=1.0\n");
  fprintf(fp,"BFGS_MaxUpdates=5\n");

  fflush(fp);
  iwave_fclose(fp);

  fprintf(stdout,"return value from creatargs = %d\n",ps_createargs(par,xargc-1,xargv+1));

  fprintf(stdout,"created PARARRAY:\n");
  ps_printall(*par,stdout);

  free(xargv[1]);
  free(xargv[2]);
  free(xargv[3]);
  free(xargv);
  ps_delete(&par);
  iwave_fdestroy();
}
