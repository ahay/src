#include <parser.h>

int main(int argc, char ** argv) {

  FILE * fp = NULL;
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

  // create test par file
  char * tname = (char *)malloc(128*sizeof(char));
  memset(tname,'\0',128);
  strcpy(tname,"test/partest/test.par");

  if (!(fp = iwave_fopen(&tname,"w",NULL,stderr))) {
    fprintf(stderr,"PANIC!\n");
    exit(1);
  }
  fprintf(fp,"n1=721 d1=2.500000e+00 o1=0.000000e+00\n");
  fprintf(fp,"n2=3121 d2=2.500000e+00 o2=0.000000e+00\n");
  fprintf(fp,"n3=1 d3=1.000000e+00 o3=3.300000e+03\n");
  fprintf(fp,"data_type=velocity\n");
  fprintf(fp,"data_format=native_float\n");
  fprintf(fp,"in=../../model/data/vp2d_2.5m.rsf@\n");
  fprintf(fp,"1.4-svn	sfdd	demo/model/data:	sergey@Sergeys-MacBook-Air.local	Thu Jun  7 02:38:15 2012\n");
  fprintf(fp,"\n");
  fprintf(fp,"	data_format=\"xdr_float\"\n");
  fprintf(fp,"	esize=4\n");
  fprintf(fp,"	in=\"stdout\"\n");
  fprintf(fp,"	in=\"stdin\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"1.4-svn	sfdd	trip/iwave/demo1:	wsymes@william-symess-macbook-pro-15.local	Sun Jul  1 10:22:23 2012\n");
  fprintf(fp,"\n");
  fprintf(fp,"	data_format=\"native_float\"\n");
  fprintf(fp,"	esize=4\n");
  fprintf(fp,"	in=\"stdout\"\n");
  fprintf(fp,"	in=\"stdin\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"1.4-svn	sfwindow	trip/iwave/demo1:	wsymes@william-symess-macbook-pro-15.local	Sun Jul  1 10:22:23 2012\n");
  fprintf(fp,"\n");
  fprintf(fp,"	d2=20.0000012314156\n");
  fprintf(fp,"	n2=391\n");
  fprintf(fp,"	o1=0\n");
  fprintf(fp,"	o2=0\n");
  fprintf(fp,"	data_format=\"native_float\"\n");
  fprintf(fp,"	esize=4\n");
  fprintf(fp,"	in=\"stdout\"\n");
  fprintf(fp,"	d1=20\n");
  fprintf(fp,"	n1=91\n");
  fprintf(fp,"	in=\"stdin\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"1.4-svn	sfput	trip/iwave/demo1:	wsymes@william-symess-macbook-pro-15.local	Sun Jul  1 10:22:23 2012\n");
  fprintf(fp,"\n");
  fprintf(fp,"	data_format=\"native_float\"\n");
  fprintf(fp,"	unit=km/s\n");
  fprintf(fp,"	label=Velocity\n");
  fprintf(fp,"	in=\"/var/tmp/trip/iwave/demo1/vp3.rsf@\"\n");
  fprintf(fp,"	unit1=m\n");
  fprintf(fp,"	label1=\n");
  fprintf(fp,"	unit2=m\n");
  fprintf(fp,"	label2=Distance\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  
  fflush(fp);
  iwave_fclose(fp);

  if (err=ps_createfile(par,"test/partest/test.par")) {
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
