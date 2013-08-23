#include <parser.h>

int main(int argc, char** argv) 
{
  PARARRAY * par = NULL;

  if (argc != 2) {
    fprintf(stderr,"Usage: ps_strip.x <file name>, writes to stdout\n");	
    exit(1);
  }

  par = ps_new();
  if (!par) {
    fprintf(stderr,"Error: ps_strip.x\n");
    fprintf(stderr,"failed to allocate new PARARRAY\n");
    exit(1);
  }
  if (ps_createfile(par,argv[1])) {
    fprintf(stderr,"Error: ps_strip.x\n");
    fprintf(stderr,"failed to initialize PARARRAY from file %s\n",argv[1]);
    exit(1);
  }
  ps_printall(*par,stdout);
  ps_delete(&par);
}
    
