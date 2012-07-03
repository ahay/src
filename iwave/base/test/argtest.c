#include <parser.h>

int main(int argc, char ** argv) {

  PARARRAY * par = ps_new();

  // create fake args
  int xargc = 4;
  char ** xargv = (char **)malloc(xargc*sizeof(char *));
  xargv[0]=NULL;
  xargv[1]=(char *)malloc(strlen("DispFlag=2")+1);
  strcpy(xargv[1],"DispFlag=0");
  xargv[2]=(char *)malloc(strlen("ResTol=0.001")+1);
  strcpy(xargv[2],"ResTol=0.001");
  xargv[3]=(char *)malloc(strlen("par=test/umin.par")+1);
  strcpy(xargv[3],"par=test/umin.par");

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
