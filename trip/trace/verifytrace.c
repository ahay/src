/* Trace comparison. */

/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#include <rsf_su.h>

#include "verifytrace.h"
#include "lci.h"

/*#define DBG*/

int verifytrace(segy * trial, /* input trial trace */ 
		segy * comp,  /* comparison trace */
		int dim,      /* dimension of wave propagation */
		float tol,    /* tolerance for success */
		float rc,     /* nominal reflection coefficient */
		float cref,   /* nominal velocity */
		float amp,    /* nominal amplitude at reference distance */
                float rdist,  /* nominal reference distance for amplitue normalization */
		float * e)   /* scaled max error */
/*< trace comparison. computes nominal amplitude curve based on
  approximate reflection amplitude as function of source and receiver
  positions and time, returns relative max norm, normalized by this
  nominal amplitude series. The nominal reflection coefficient used
  to compute this reference amplitude is the fourth argument, rc.
  *
  The direct wave amplitude, is amp*rdist/dist, where rdist is a
  reference distance - thus amp is meant to be the amplitude of the
  direct wave at the reference distance. The direct wave amplitude,
  at a distance equal to the two-way time multiplied by reference
  velocity, scaled by the refl coeff, is the proxy for reflected wave
  reference amplitude.
  *
  To accomodate possibly different time steps and numbers of samples,
  use local cubic interpolation.
  >*/
{
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
