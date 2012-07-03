#include <parser.h>

int main(int argc, char ** argv) {

  int err=0;
  PARARRAY * par = ps_new();
  char * s = NULL;
  char * t = NULL;
  char * dff = NULL;
  char * inf = NULL;
  char * dfl = NULL;
  char * inl = NULL;
  char * dffa = NULL;
  char * dfla = NULL;
  long l = 0L;
  short h = 0;
  int i = 0;
  float x = 0.0f;
  double xx = 0.0;

  if (err=ps_createfile(par,"test/test.par")) {
    fprintf(stderr,"Error from ps_createfile: err=%d\n",err);
    ps_delete(&par);
    exit(1);
  }

  printf("DUMP of PARARRAY, before:\n");
  printf("------------------------------------------------------\n");
    
  ps_printall(*par,stdout);
  printf("------------------------------------------------------\n\n");
  if (err=ps_ffcstring(*par,"n1",&s)) {
    fprintf(stdout,"did not find n1 (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of n1 (cstring) = %s\n",s);
  }
  if (err=ps_ffcstring(*par,"fruitcake",&t)) {
    fprintf(stdout,"did not find fruitcake (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of fruitcake (cstring) = %s\n",t);
  }

  if (err=ps_fflong(*par,"n1",&l)) {
    fprintf(stdout,"did not find n1 (long), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of n1 (long) = %ld\n",l);
  }

  if (err=ps_ffcstring(*par,"data_format",&dff)) {
    fprintf(stdout,"did not find data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of data_format (cstring) = %s\n",dff);
  }

  if (err=ps_flcstring(*par,"data_format",&dfl)) {
    fprintf(stdout,"did not find data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"last value of data_format (cstring) = %s\n",dfl);
  }

  if (err=ps_fflong(*par,"data_format",&l)) {
    fprintf(stdout,"did not find data_format (long), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of data_format (long) = %ld\n",l);
  }

  if (err=ps_ffcstring(*par,"in",&inf)) {
    fprintf(stdout,"did not find in (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of in (cstring) = %s\n",inf);
  }

  if (err=ps_flcstring(*par,"in",&inl)) {
    fprintf(stdout,"did not find in (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"last value of in (cstring) = %s\n",inl);
  }

  if (err=ps_fflong(*par,"d1",&l)) {
    fprintf(stdout,"did not find d1 (long), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of d1 (long) = %ld\n",l);
  }

  if (err=ps_ffshort(*par,"n2",&h)) {
    fprintf(stdout,"did not find n2 (short), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of n2 (short) = %hd\n",h);
  }

  if (err=ps_ffint(*par,"n2",&i)) {
    fprintf(stdout,"did not find n2 (int), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of n2 (int) = %d\n",h);
  }

  if (err=ps_fffloat(*par,"n2",&x)) {
    fprintf(stdout,"did not find n2 (float), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of n2 (float) = %f\n",x);
  }
  if (err=ps_fffloat(*par,"d2",&x)) {
    fprintf(stdout,"did not find d2 (float), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of d2 (float) = %10.6f\n",x);
  }  

  if (err=ps_ffdouble(*par,"d2",&xx)) {
    fprintf(stdout,"did not find d2 (double), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of d2 (double) = %20.14lf\n",xx);
  }  

  if (err=ps_fldouble(*par,"d2",&xx)) {
    fprintf(stdout,"did not find d2 (double), err=%d\n",err);
  }
  else {
    fprintf(stdout,"last value of d2 (double) = %20.14lf\n",xx);
  }  

  // change values of data_format
  if (err=ps_sfcstring(*par,"data_format","xdr_float")) {
    fprintf(stdout,"did not set data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"set first value of data_format (cstring) = %s\n","xdr_float");
  }

  if (err=ps_slcstring(*par,"data_format","xdr_float")) {
    fprintf(stdout,"did not set data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"set last value of data_format (cstring) = %s\n","xdr_float");
  }
  
  if (err=ps_ffcstring(*par,"data_format",&dffa)) {
    fprintf(stdout,"did not find data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"first value of data_format (cstring) = %s\n",dffa);
  }

  if (err=ps_flcstring(*par,"data_format",&dfla)) {
    fprintf(stdout,"did not find data_format (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"last value of data_format (cstring) = %s\n",dfla);
  }

  // set a couple more things just for the heck of it
  if (err=ps_slcstring(*par,"unit2","m")) {
    fprintf(stdout,"did not set unit2 (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"set last value of unit2 (cstring) = %s\n","m");
  }
  if (err=ps_slcstring(*par,"label2","Depth")) {
    fprintf(stdout,"did not set label2 (cstring), err=%d\n",err);
  }
  else {
    fprintf(stdout,"set last value of label2 (cstring) = %s\n","Depth");
  }

  printf("DUMP of PARARRAY, after:\n");
  printf("------------------------------------------------------\n");
    
  ps_printall(*par,stdout);
  printf("------------------------------------------------------\n\n");

  if (s) userfree_(s);
  if (t) userfree_(t);
  if (dff) userfree_(dff);
  if (dfl) userfree_(dfl);
  if (inf) userfree_(inf);
  if (inl) userfree_(inl);
  if (dffa) userfree_(dffa);
  if (dfla) userfree_(dfla);
  ps_delete(&par);

  iwave_fdestroy();
}
