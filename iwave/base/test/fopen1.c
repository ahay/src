#include "iwave_fopen.h"

int main(int argc, char ** argv) {

  FILE * fp1 = NULL;
  FILE * fp2 = NULL;
  FILE * fp3 = NULL;
  FILE * fp4 = NULL;
  FILE * fp5 = NULL;
  FILE * fp6 = NULL;
  FILE * fp7 = NULL;
  FILE * fp8 = NULL;

  system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > h1.su");
  system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > f1.su");
  system("sunull nt=201 ntr=21 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > f2.su");

  char * buf = (char *)malloc(128*sizeof(char));
  strcpy(buf,"h1.su");
  fp1=iwave_fopen(&buf,"r",NULL,stdout);
  strcpy(buf,"f1.su");
  fp2=iwave_fopen(&buf,"r","h1.su",stdout);
  strcpy(buf,"n1.su");
  fp3=iwave_fopen(&buf,"w+","h1.su",stdout);
  char * blank=NULL;
  fp4=iwave_fopen(&blank,"w+","h1.su",stdout);
  strcpy(buf,"f2.su");
  fp5=iwave_fopen(&buf,"r",NULL,stdout);
  strcpy(buf,"h1.su");
  fp6=iwave_fopen(&buf,"r+","h1.su",stdout);

  fprintf(stdout,"AFTER INITIAL OPENS\n");

  iwave_fprintall(stdout);

  free(blank);
  blank=NULL;
  fp6=iwave_fopen(&blank,"w+","h1.su",stdout);

  iwave_fclose(fp4);
  fprintf(stdout,"AFTER CLOSE OF TEMP FILE\n");
  iwave_fprintall(stdout);

  free(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+","h1.su",stdout);
  free(blank);
  blank=NULL;
  fp4=iwave_fopen(&blank,"w+","f2.su",stdout);

  fprintf(stdout,"OPEN OF TWO MORE\n");
  fprintf(stdout,"ONE ON SAME PROTOTYPE, ONE ON DIFFERENT\n");

  iwave_fprintall(stdout);

  fprintf(stdout,"ATTEMPT ERROR CONDITIONS:\n");

  iwave_fclose(fp7);
  free(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"r","f2.su",stdout);

  fprintf(stdout,"attempt temp open with r access, got FILE * = %p\n",(void *)fp7);

  iwave_fclose(fp7);
  free(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+",NULL,stdout);
  
  fprintf(stdout,"attempt to open tmp file without proto, got FILE * = %p\n",(void *)fp7);

  fprintf(stdout,"AFTER ERROR CONDITIONS:\n");

  iwave_fprintall(stdout);
 
  fprintf(stdout,"REQUEST FILE PTR FOR TMP FILE TEMPLATED ON h1\n");

  free(blank);
  blank=NULL;
  fp8=iwave_fopen(&blank,"w+","h1.su",stdout);

  free(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+","h1.su",stdout);

  iwave_fclose(fp4);
  
  free(blank);
  blank=NULL;
  fp4=iwave_fopen(&blank,"w+","h1.su",stdout);

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

  free(blank);
  blank=NULL;
  fp7=iwave_fopen(&blank,"w+","h1.su",stdout);

  free(blank);
  blank=NULL;
  fp6=iwave_fopen(&blank,"w+","f2.su",stdout);

  fprintf(stdout,"OPEN TEMP ON FP7 PROTO H1, FP6 PROTO F2\n");

  iwave_fprintall(stdout);
 
  fprintf(stdout,"TRASH IT ALL\n");
  
  iwave_fdestroy();
  iwave_fprintall(stdout);

  return(0);
}  

  

