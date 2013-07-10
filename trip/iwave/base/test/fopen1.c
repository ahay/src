#include "iwave_fopen.h"

#define FILE1 "file1"
#define FILE2 "file2"
#define FILE3 "file3"

int main(int argc, char ** argv) {

    /* generate 3 data files */
  system("echo \"sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100\" > file1");
  system("echo \"sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100\" > file2");
  system("echo \"sunull nt=201 ntr=21 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100\" > file3");

  FILE * fp1 = NULL;
  FILE * fp2 = NULL;
  FILE * fp3 = NULL;
  FILE * fp4 = NULL;
  FILE * fp5 = NULL;
  FILE * fp6 = NULL;
  FILE * fp7 = NULL;
  FILE * fp8 = NULL;
  FILE * fp9 = NULL;

  char * buf = NULL;
  char * tn1 = NULL;
  char * tn2 = NULL;
  char * tn3 = NULL;
  char * tn4 = NULL;
  char * tn5 = NULL;
  char * tn6 = NULL;
  char * blank = NULL;

  fprintf(stdout,"OPEN:\n");
  fprintf(stdout,"old file no proto: %s mode=r proto=NULL\n",FILE1);
  buf = (char *)usermalloc_(128*sizeof(char));
  strcpy(buf,FILE1);
  fp1=iwave_fopen(&buf,"r",NULL,stdout);
  strcpy(buf,FILE2);
  fprintf(stdout,"old file with proto: %s mode=r proto=%s\n",FILE2,FILE1);
  fp2=iwave_fopen(&buf,"r",FILE1,stdout);
  strcpy(buf,"n1.su");
  fprintf(stdout,"new file with proto: n1.su mode=w+ proto=%s\n",FILE1);
  fp3=iwave_fopen(&buf,"w+",FILE1,stdout);
  fprintf(stdout,"tmp file with proto: ");
  fp4=iwave_fopen(&tn1,"w+",FILE1,stdout);
  fprintf(stdout,"%s mode=w+ proto=%s\n",tn1,FILE1);
  fprintf(stdout,"old file without proto: %s mode=r proto=NULL\n",FILE3);
  strcpy(buf,FILE3);
  fp5=iwave_fopen(&buf,"r",NULL,stdout);
  fprintf(stdout,"old file, reopened with r rather than r+ mode: %s mode=r+ proto=NULL\n",FILE1);
  strcpy(buf,FILE1);
  fp6=iwave_fopen(&buf,"r+",NULL,stdout);
  if (fp6==fp1) fprintf(stdout," -- returned same pointer\n");

  fprintf(stdout,"AFTER INITIAL OPENS\n");

  iwave_fprintall(stdout);

  fp6=iwave_fopen(&tn2,"w+",FILE1,stdout);
  fprintf(stdout,"tmp file with proto: %s mode=w+ proto=%s\n",tn2,FILE1);

  iwave_fclose(fp4);
  fprintf(stdout,"close %s\n",tn1);
  fprintf(stdout,"AFTER OPEN OF NEW TEMP FILE PROTO, CLOSE OF FIRST TEMP FILE\n");
  iwave_fprintall(stdout);

  fp7=iwave_fopen(&tn3,"w+",FILE1,stdout);
  fprintf(stdout,"\n\ntmp file with proto: %s mode=w+ proto=%s\n",tn3,FILE1);
  fprintf(stdout,"should return %s since same proto and latter is closed (not in use)\n",tn1);
  if (fp7 == fp4) fprintf(stdout," -- returned same pointer\n");

  fprintf(stdout,"attempt to open tmp file %s on another  FILE*, should retrieve pointer\n",tn3);
  fprintf(stdout,"use mode=r to test that it's treated as open anyway\n");
  fp9=iwave_fopen(&tn3,"r",FILE1,stdout);
  if (fp9==fp7) fprintf(stdout," -- returned same pointer\n");

  iwave_fprintall(stdout);

  fprintf(stdout,"\n\n now close duplicate pointer\n");
  iwave_fclose(fp9);
  iwave_fprintall(stdout);

  fprintf(stdout,"\n\n then close original pointer\n");
  iwave_fclose(fp7);
  iwave_fprintall(stdout);

  fp7=iwave_fopen(&tn4,"w+",FILE3,stdout);
  fprintf(stdout,"\,\,tmp file with proto: %s mode=w+ proto=%s\n",tn4,FILE3);
  fprintf(stdout,"should NOT return %s since protos are different\n",tn1);

  fprintf(stdout,"AFTER TEMP FILE MANIPULATIONS:\n");

  iwave_fprintall(stdout);

  fprintf(stdout,"ATTEMPT ERROR CONDITIONS:\n");
  fprintf(stdout,"close %s\n",tn4);
  iwave_fclose(fp7);
  fp7=NULL;
  fp7=iwave_fopen(&blank,"r",FILE3,stdout);

  fprintf(stdout,"attempt temp open with r access, got FILE * = %p\n",(void *)fp7);

  fp7=NULL;
  fp7=iwave_fopen(&blank,"w+",NULL,stdout);
  
  fprintf(stdout,"attempt to open tmp file without proto, got FILE * = %p\n",(void *)fp7);

  fprintf(stdout,"AFTER ERROR CONDITIONS:\n");

  iwave_fprintall(stdout);
 
  fprintf(stdout,"REQUEST 2 FILE PTRS FOR TMP FILES TEMPLATED ON %s\n",FILE1);

  fp8=iwave_fopen(&tn5,"w+",FILE1,stdout);
  fp9=iwave_fopen(&tn6,"w+",FILE1,stdout);

  iwave_fprintall(stdout);

  fprintf(stdout,"CLOSE %s, REOPEN\n",tn1); 
  iwave_fclose(fp4);
  userfree_(tn1);
  tn1=NULL;
  fp4=iwave_fopen(&tn1,"w+",FILE1,stdout);

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
  iwave_fclose(fp9);

  fprintf(stdout,"EVERYTHING IWAVE_FCLOSED\n");

  iwave_fprintall(stdout);

  fprintf(stdout,"TRASH IT ALL\n");

  userfree_(buf);
  userfree_(tn1);
  userfree_(tn2);
  userfree_(tn3);
  userfree_(tn4);
  userfree_(tn5);
  userfree_(tn6);

  iwave_fdestroy();
  iwave_fprintall(stdout);

  return(0);
}  

  

