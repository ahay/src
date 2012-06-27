#include <traceio.h>
#include <parser.h>

#define FILE1 "test/data0.su"
#define FILE2 "test/data3.su"

char ** xargv;

int main(int argc, char ** argv) {

  segy trial;
  segy comp;
  FILE * fptrial=NULL;
  FILE * fpcomp=NULL;
  ireal emax = REAL_ZERO;
  ireal e;
  int i;
  int ntr = 0;

  char * trialname = (char *)usermalloc_(sizeof(char)*(strlen(FILE1)+1));
  char * compname = (char *)usermalloc_(sizeof(char)*(strlen(FILE2)+1));

  if (!trialname || !compname) {
    fprintf(stderr,"Error:trcomp\n");
    fprintf(stderr,"failed to allocate filename arrays\n");
    exit(1);
  }
  strcpy(trialname,FILE1);
  strcpy(compname,FILE2);

  if (!(fptrial=iwave_fopen(&trialname,"r",NULL,stderr)) ||
      !(fpcomp =iwave_fopen(&compname,"r",NULL,stderr))) {
    fprintf(stderr,"Error: trcomp\n");
    fprintf(stderr,"failed to open files\n");
    exit(1);    
  }

  userfree_(trialname);
  userfree_(compname);

  while (fgettr(fptrial,&trial) && fgettr(fpcomp,&comp)) {
    if (trial.ns != comp.ns) {
      fprintf(stderr,"Error: trcomp\n");
      fprintf(stderr,"traces have unequal length\n");
      exit(1);
    }
    e=REAL_ZERO;
    for (i=0;i<trial.ns;i++) 
      e=iwave_max(e,iwave_abs(trial.data[i]-comp.data[i]));

    printf("trace %d: max abs diff = %e\n", 
	   ntr,e);
    emax=fmax(e,emax);
    ntr++;
  }
  printf("max error over all traces = %e\n",emax);
  iwave_fclose(fptrial);
  iwave_fclose(fpcomp);

  iwave_fdestroy();
}
