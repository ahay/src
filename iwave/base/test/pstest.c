#include <keyval.h>

int main(int argc, char ** argv) {


  PSLINK * par = pslink_new();
  WORD * key = word_new();
  WORD * val = word_new();
  char * teststr = (char *)malloc(1000*sizeof(char));

  strcpy(teststr,
	 "a test of the pair reader  \
          \n just to see if it works \na\n=  \
          b k=wwww0w \" this is a\n quote \" \
          this=fun this  =     \n \"no fun\" \
          ha ha ha mary had a little=lamb===");
  /*
  do {
    kv_read(par->pair,&teststr);
    stop = kv_check(*(par->pair));
    if (!stop) {
      fprintf(stderr,"pstest: added ");
      kv_fprint(*(par->pair),stderr);
      par->next = pslink_new();
      par->next->prev = par;
      par       = par->next;
    }
    else {
      par = par->prev;
      pslink_delete(&(par->next));
    }
  } while (!stop);
  */

  if (pslink_read(&par,&teststr)) {
    fprintf(stderr,"Error: pstest - failed to read data\n");
    exit(1);
  }

  fprintf(stderr,"\nback to root\n");
  //  while (par->prev) par=par->prev;
  if (pslink_front(&par)) {
    fprintf(stderr,"Error: pstest - failed to move to front\n");
    exit(1);
  }

  fprintf(stderr,"\nprint in normal order\n");
  while (par->next) {
    kv_fprint(*(par->pair),stderr);
    par = par->next;
  }
  kv_fprint(*(par->pair),stderr);
  
  fprintf(stderr,"\nprint in reverse order\n");
  while (par->prev) {
    kv_fprint(*(par->pair),stderr);
    par = par->prev;
  }
  kv_fprint(*(par->pair),stderr);

  fprintf(stderr,"\nfind first instance of key=\"this\", print value\n");
  word_assign(key,"this",4);
  pslink_findfirst(par,*key,val);
  fprintf(stderr,"\nkey = this val = ");
  fprintf(stderr,"%s\n",val->str);

  fprintf(stderr,"\nreset value for existing key\n");
  word_assign(val,"\"totally fun\"",strlen("\"totally fun\""));
  pslink_setfirst(&par,*key,*val);

  fprintf(stderr,"\nset value for new key\n");
  word_assign(key,"that",4);
  word_assign(val,"\"not much fun at all\"",22);
  pslink_setfirst(&par,*key,*val);

  fprintf(stderr,"\nto the front\n");
  pslink_front(&par);
  fprintf(stderr,"\nprint again in normal order\n");
  while (par->next) {
    kv_fprint(*(par->pair),stderr);
    par = par->next;
  }
  kv_fprint(*(par->pair),stderr);

}
