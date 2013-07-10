#include <parser.h>

int main(int argc, char ** argv) {

  WORD * w = word_new();

  int wc=0;

  char * teststr = (char *)malloc(1000*sizeof(char));
  char * finger = teststr;
  strcpy(teststr,
	 "a test of the word reader  \
          \n just to see if it works \na\n=  \
          b k=wwww0w \" this is a\n quote \" \
          this=fun this  =     \n \"no fun\" \
          ha ha ha mary had a little=lamb=== \
           a=c=b=d");

  fprintf(stdout,"%s\n",teststr);

  fprintf(stdout,"\n --- here are the words:\n");
  
  while (!word_read(w,&finger) && w->str) {
    fprintf(stdout,"word %d = %s\n",wc,w->str);
    wc++;
  }

  word_delete(&w);

  free(teststr);
}
