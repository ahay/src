#include "iwave_fopen.h"

#define FILE1 "test/h1.su"
#define FILE2 "test/f1.su"
#define FILE3 "test/f2.su"

int main(int argc, char ** argv) {

  FILE * fp1 = NULL;
  FILE * fp2 = NULL;
  FILE * fp3 = NULL;
  FILE * fp4 = NULL;
  FILE * fp5 = NULL;
  FILE * fp6 = NULL;
  FILE * fp7 = NULL;
  FILE * fp8 = NULL;

  char * buf = (char *)usermalloc_(128*sizeof(char));
  strcpy(buf,FILE1);
  fp1=iwave_fopen(&buf,"r",NULL,stdout);
  strcpy(buf,FILE2);
  fp2=iwave_fopen(&buf,"r",FILE1,stdout);
  strcpy(buf,"test/fopen1/n1.su");
  fp3=iwave_fopen(&buf,"w+",FILE1,stdout);
  char * blank=NULL;
  fp4=iwave_fopen(&blank,"w+",FILE1,stdout);
  strcpy(buf,FILE3);
  fp5=iwave_fopen(&buf,"r",NULL,stdout);
  strcpy(buf,FILE1);
  fp6=iwave_fopen(&buf,"r+",FILE1,stdout);

  fprintf(stdout,"AFTER INITIAL OPENS\n");

  iwave_fprintall(stdout);

  userfree_(blank);
  blank=NULL;
  fp6=iwave_fopen(&blank,"w+",FILE1,stdout);

  iwave_fclose(fp4);
  fprintf(stdout,"AFTER CLOSE OF TEMP FILE\n");
  iwave_fprintall(stdout);

  userfree_(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+",FILE1,stdout);
  userfree_(blank);
  blank=NULL;
  fp4=iwave_fopen(&blank,"w+",FILE3,stdout);

  fprintf(stdout,"OPEN OF TWO MORE\n");
  fprintf(stdout,"ONE ON SAME PROTOTYPE, ONE ON DIFFERENT\n");

  iwave_fprintall(stdout);

  fprintf(stdout,"ATTEMPT ERROR CONDITIONS:\n");

  iwave_fclose(fp7);
  userfree_(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"r",FILE3,stdout);

  fprintf(stdout,"attempt temp open with r access, got FILE * = %p\n",(void *)fp7);

  iwave_fclose(fp7);
  userfree_(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+",NULL,stdout);
  
  fprintf(stdout,"attempt to open tmp file without proto, got FILE * = %p\n",(void *)fp7);

  fprintf(stdout,"AFTER ERROR CONDITIONS:\n");

  iwave_fprintall(stdout);
 
  fprintf(stdout,"REQUEST FILE PTR FOR TMP FILE TEMPLATED ON h1\n");

  userfree_(blank);
  blank=NULL;
  fp8=iwave_fopen(&blank,"w+",FILE1,stdout);

  userfree_(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+",FILE1,stdout);

  iwave_fclose(fp4);
  
  userfree_(blank);
  blank=NULL;
  fp4=iwave_fopen(&blank,"w+",FILE1,stdout);

  fprintf(stdout,"AFTER MORE OPENS: fp4=%p fp7=%p fp8=%p\n",(void *)fp4,(void *)fp7,(void *)fp8);

  iwave_fprintall(stdout);

  iwave_fclose(fp1);
  iwave_fclose(fp2);
  iwave_fclose(fp3);
  iwave_fclose(fp4);
  iwave_fclose(fp5);
  iwave_fclose(fp6);
  iwave_fclose(fp7);
  iwave_fclose(fp8);

  fprintf(stdout,"EVERYTHING IWAVE_FCLOSED\n");

  iwave_fprintall(stdout);

  userfree_(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+",FILE1,stdout);

  userfree_(blank);
  blank=NULL;
  fp6=iwave_fopen(&blank,"w+",FILE3,stdout);

  fprintf(stdout,"OPEN TEMP ON FP7 PROTO H1, FP6 PROTO F2\n");

  iwave_fprintall(stdout);
 
  fprintf(stdout,"TRASH IT ALL\n");

  userfree_(blank);
  userfree_(buf);
  iwave_fdestroy();
  iwave_fprintall(stdout);

  return(0);
}  

  

