#include <keyval.h>

int main(int argc, char ** argv) {


  PSLINK * par = pslink_new();
  WORD * key = word_new();
  WORD * val = word_new();
  char * teststr = (char *)malloc(1000*sizeof(char));
  char * save = teststr;

  strcpy(teststr,
	 "a test of the pair reader  \
          \n just to see if it works \na\n=  \
          b k=wwww0w \" this is a\n quote \" \
          this=fun this  =     \n \"no fun\" \
          \"quoted = equation\"    \
          ha ha ha mary had a little=lamb===\
          \n    a=c=b=d");

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

  fprintf(stderr,"\nreset value for existing key - first\n");
  word_assign(val,"\"totally fun\"",strlen("\"totally fun\""));
  pslink_setfirst(&par,*key,*val);

  fprintf(stderr,"\nreset value for existing key - last\n");
  word_assign(val,"\"bummer\"",strlen("\"bummer\""));
  pslink_setlast(&par,*key,*val);

  fprintf(stderr,"\nset value for new key - first\n");
  word_assign(key,"that",4);
  word_assign(val,"\"not much fun at all\"",22);
  pslink_setfirst(&par,*key,*val);

  fprintf(stderr,"\nset value for new key - last\n");
  word_assign(key,"\"the other\"",13);
  word_assign(val,"\"quite a lot of fun actually\"",29);
  pslink_setlast(&par,*key,*val);

  fprintf(stderr,"\nto the front\n");
  pslink_front(&par);
  fprintf(stderr,"\nprint again in normal order\n");
  while (par->next) {
    kv_fprint(*(par->pair),stderr);
    par = par->next;
  }
  kv_fprint(*(par->pair),stderr);

  fprintf(stderr,"\nfind first instance of key \"a\"\n");
  word_assign(key,"a",4);
  pslink_findfirst(par,*key,val);
  fprintf(stderr,"\nkey = %s val = %s\n",key->str,val->str);
  fprintf(stderr,"\nfind last instance of key \"a\"\n");
  word_assign(key,"a",4);
  pslink_findlast(par,*key,val);
  fprintf(stderr,"\nkey = %s val = %s\n",key->str,val->str);

  fprintf(stderr,"\nre-read string, should add same pairs all over again at end\n");
  teststr = save;
  pslink_read(&par,&teststr);
  pslink_front(&par);
  fprintf(stderr,"\nprint again in normal order\n");
  while (par->next) {
    kv_fprint(*(par->pair),stderr);
    par = par->next;
  }
  kv_fprint(*(par->pair),stderr);

  free(save);
  word_delete(&val);
  word_delete(&key);
  pslink_delete(&par);
}
