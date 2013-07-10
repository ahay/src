#include <parser.h>

int main(int argc, char** argv) {

  if (argc != 3) {
    fprintf(stderr,"Usage: ps_getstring.x <file name> <key>, writes val to stdout\n");	
    fprintf(stderr,"");
    exit(1);
  }
  PARARRAY * par = NULL;
  par = ps_new();
  if (!par) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to allocate new PARARRAY\n");
    fprintf(stderr,"");
    exit(1);
  }
  if (ps_createfile(par,argv[1])) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to initialize PARARRAY from file %s\n",argv[1]);
    fprintf(stderr,"");
    exit(1);
  }
  char * val = NULL;
  ps_flcstring(*par,argv[2],&val);
  if (!val) {
    fprintf(stderr,"Error: ps_getstring.x\n");
    fprintf(stderr,"failed to extract key %s PARARRAY from file %s\n",argv[2],argv[1]);
    fprintf(stdout,"");
  }
  else {
    fprintf(stdout,"%s",val);
    userfree_(val);
  }
  ps_delete(&par);
}
    
