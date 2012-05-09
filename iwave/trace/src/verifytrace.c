#include "verifytrace.h"

/*#define DBG*/

int verifytrace(segy * trial, 
		segy * comp,  
		int dim,
		float tol,    
		float rc,     
		float cref,  
		float amp,
		float rdist,
		float * e) {  

  int i;
  float fac;
  float le;
  float lemax;
  float tdt=0.001*(float)(trial->dt);
  float cdt=0.001*(float)(comp->dt);
  int tnt=trial->ns;
  int cnt=comp->ns;
  float tot=trial->delrt;
  float cot=comp->delrt;
  float tmax=0.0;
  float cmax=0.0;
  
  float *ctimes = (float*)malloc(cnt*sizeof(float));
  float *tint   = (float*)malloc(cnt*sizeof(float));

  double oops;
  double duh;

#ifdef DBG
  FILE * fp;
  segy t;
#endif

  for (i=0;i<cnt;i++) ctimes[i]=cot+i*cdt;

  lci(cnt,tnt,tdt,tot,ctimes,trial->data,tint);

#ifdef DBG
  t.tracl=1;
  t.ns=cnt;
  t.dt=(int)(1000*cdt);
  for (i=0;i<cnt;i++) {
    t.data[i]=tint[i];
  }
  fp=fopen("testver.su","w");
  fputtr(fp,&t);
  fclose(fp);
#endif
  
  *e=0.0;


  lemax=0.0;
  if (dim==3) {
    fac=rc*amp*rdist/cref;
    for (i=0;i<cnt;i++) {
      oops=tint[i]-comp->data[i];
      duh=fabs(oops);
      le=duh;
      /* cannot use in c90
      lemax=fmaxf(le,lemax);
      */
      oops=lemax;
      lemax=fmax(duh,oops);
      oops=*e;
      duh=le*(cot+i*cdt)/fac;
      /*
      *e=fmaxf(*e,
	      le*(cot+i*cdt)/fac);
      */
      *e=fmax(oops,duh);
      /*
      cmax=fmaxf(fabsf(comp->data[i]),cmax);
      */
      oops=cmax;
      duh=comp->data[i];
      duh=fmax(oops,fabs(duh));
      cmax=duh;
      /*
      tmax=fmaxf(fabsf(tint[i]),tmax);
      */
      oops=tmax;
      duh=tint[i];
      duh=fmax(fabs(duh),oops);
      tmax=duh;
#ifdef DBG
      printf("it relerr,cumrelerr,comperr,tgtabs,compabs,abserr,cumabserr\n");
      printf("%d %e %e %e %e %e %e\n",i,fabsf((tint[i]-comp->data[i])*(cot+i*cdt)/fac),*e,fac/(cot+i*cdt),fabsf(comp->data[i]),le,lemax);
#endif
    }
  }
  else if (dim==2) {
    fac=rc*amp*sqrt(rdist/cref);
    for (i=0;i<cnt;i++) {
      oops=tint[i]-comp->data[i];
      duh=fabs(oops);
      le=duh;
      /* cannot use in c90
      lemax=fmaxf(le,lemax);
      */
      oops=lemax;
      lemax=fmax(duh,oops);
      oops=*e;
      duh=cot+i*cdt;
      duh=sqrt(fabs(duh));
      duh=duh*le/fac;
      /*
      *e=fmaxf(*e,
	      le*sqrtf(fabsf(cot+i*cdt))/fac);
      */
      *e=fmax(oops,duh);
      /*
      cmax=fmaxf(fabsf(comp->data[i]),cmax);
      */
      oops=cmax;
      duh=comp->data[i];
      duh=fmax(oops,fabs(duh));
      cmax=duh;
      /*
      tmax=fmaxf(fabsf(tint[i]),tmax);
      */
      oops=tmax;
      duh=tint[i];
      duh=fmax(fabs(duh),oops);
      tmax=duh;
#ifdef DBG
      printf("it relerr,cumrelerr,comperr,tgtabs,compabs,abserr,cumabserr\n");
      printf("%d %e %e %e %e %e %e\n",i,fabsf((tint[i]-comp->data[i])*(cot+i*cdt)/fac),*e,fac/(cot+i*cdt),fabsf(comp->data[i]),le,lemax);
#endif
    }
  }
  else {
    fprintf(stderr,"Error: verify - dim must be 2 or 3\n");
    exit(1);
  }
    
#ifdef DBG
  printf("max of test trace = %e\n",tmax);
  printf("max of comp trace = %e\n",cmax);
#endif
  free(ctimes);
  free(tint);
  if (*e<tol) return 1;
  return 0;
}

  /*
  Value val;

  float tsz,tsx,tsy;
  float tgz,tgx,tgy;
  float csz,csx,csy;
  float cgz,cgx,cgy;
  float scalco, scalel;
  */
  /* compute source, receiver positions */
  /*
  gethdval(trial,"scalco",&val);
  scalco=vtof(hdtype("scalco",&val));
  gethdval(trial,"gx",&val);
  tgx = vtof(hdtype("gx"),val);
  gethdval(trial,"sx",&val);
  tsx = vtof(hdtype("sx"),val);
  gethdval(trial,"gy",&val);
  tgy = vtof(hdtype("gy"),val);
  gethdval(trial,"sy",&val);
  tsy = vtof(hdtype("sy"),val);
  if (scalco > 0) { tgx *= scalco; tgy *=  scalco; }
  if (scalco < 0) { tgx /=-scalco; tgy /= -scalco; }
  if (scalco > 0) { tsx *= scalco; tsy *=  scalco; }
  if (scalco < 0) { tsx /=-scalco; tsy /= -scalco; }
  gethdval(trial,"scalel",&val);
  scalel=vtof(hdtype("scalel"),&val);
  gethdval(trial,"gelev",&val);
  tgz=-vtof(hdtype("gelev"),&val);
  gethdval(trial,"selev",&val);
  tsz=-vtof(hdtype("selev"),&val);
  if (scalel > 0) { tgz *= scalel; tsz *=  scalel; }
  if (scalel < 0) { tgz /=-scalel; tsz /= -scalel; }
  
  gethdval(comp,"scalco",&val);
  scalco=vtof(hdtype("scalco",&val));
  gethdval(comp,"gx",&val);
  cgx = vtof(hdtype("gx"),val);
  gethdval(comp,"sx",&val);
  csx = vtof(hdtype("sx"),val);
  gethdval(comp,"gy",&val);
  cgy = vtof(hdtype("gy"),val);
  gethdval(comp,"sy",&val);
  csy = vtof(hdtype("sy"),val);
  if (scalco > 0) { cgx *= scalco; cgy *=  scalco; }
  if (scalco < 0) { cgx /=-scalco; cgy /= -scalco; }
  if (scalco > 0) { csx *= scalco; csy *=  scalco; }
  if (scalco < 0) { csx /=-scalco; csy /= -scalco; }
  gethdval(comp,"scalel",&val);
  scalel=vtof(hdtype("scalel"),&val);
  gethdval(comp,"gelev",&val);
  cgz=-vtof(hdtype("gelev"),&val);
  gethdval(comp,"selev",&val);
  csz=-vtof(hdtype("selev"),&val);
  if (scalel > 0) { cgz *= scalel; csz *=  scalel; }
  if (scalel < 0) { cgz /=-scalel; csz /= -scalel; }
  */
