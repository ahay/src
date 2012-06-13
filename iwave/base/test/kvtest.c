#include <keyval.h>

int main(int argc, char ** argv) {

  KEYVAL * kv = kv_new();
  char * teststr = (char *)malloc(1000*sizeof(char));

  strcpy(teststr,
	 "a test of the pair reader  \
          \n just to see if it works \na\n=  \
          b k=wwww0w \" this is a\n quote \" \
          this=fun this  =     \n \"no fun\" \
          ha ha ha mary had a little=lamb===");
  
  do {
    kv_read(kv,&teststr);
    fprintf(stderr,"kvtest: ");
    kv_fprint(*kv,stderr);
  } while (!kv_check(*kv));

  kv_delete(&kv);
}
