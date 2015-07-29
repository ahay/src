#include <parser.h>

int main(int argc, char** argv) 
{
  PARARRAY * par = NULL;
  char * val = NULL;

  if (argc != 3) {
    fprintf(stderr,"Usage: ps_getstring.x <file name> <key>, writes val to stdout\n");	
    exit(1);
  }
  par = ps_new();
  if (!par) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to allocate new PARARRAY\n");
    exit(1);
  }
  if (ps_createfile(par,argv[1])) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to initialize PARARRAY from file %s\n",argv[1]);
    exit(1);
  }

  ps_flcstring(*par,argv[2],&val);
  if (!val) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to extract key %s PARARRAY from file %s\n",argv[2],argv[1]);
  }
  else {
    fprintf(stdout,"%s",val);
    userfree_(val);
  }
  ps_delete(&par);
  iwave_fdestroy();
}
    
