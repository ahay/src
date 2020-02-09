#include <rsf.h>

#ifndef _PREP_H

typedef struct in_para_struct in_para_struct_t;
/*^*/

typedef struct wfl_struct wfl_struct_t;
/*^*/

typedef struct acq_struct acq_struct_t;
/*^*/

typedef struct mod_struct mod_struct_t;
/*^*/

struct in_para_struct{
  bool verb;
  bool fsrf;
  bool dabc;
  bool adj;
  bool dpt;
  int nb;
};
/*^*/

struct wfl_struct{
  // pressure
  float *pc;  // current time
  float *pp;  // previous time
  float *pa;  // aux
  // velocities
  float *v1c; // component 1 current
  float *v2c; // component 2 current
  float *v1p; // component 1 previous
  float *v2p; // component 2 previous
  float *v1a; // component 1 aux
  float *v2a; // component 1 aux
  float *tap1;  // ABC coeffs
  float *tap2;  //
  // buffers to store the data
  float *rdata;
  // buffer for the wavefield snapshot
  float *bwfl;
  // dimensions
  long nabc;  // size of the absorbing boundary
  long modN1;
  long modN2;
  long simN1;
  long simN2;
  float modO1;
  float modO2;
  float simO1;
  float simO2;
  float d1;
  float d2;
  // pointers to files
  sf_file Fdata;
  sf_file Fwfl;
};
/*^*/

struct acq_struct{
  // wavelet and time parameters
  int nt;
  float dt;
  float ot;
  float *wav;
  // acquisition geometry
  long ns;    // NUMBER OF SOURCES TO BE SIMULATED SIMULTANEOUSLY
  long nr;
  float *scoord;
  float *rcoord;
  // interpolation coefficients for source injection and extraction
  float *hicksSou1;
  float *hicksSou2;
  float *hicksRcv1;
  float *hicksRcv2;
};
/*^*/

struct mod_struct{
  long n1;
  long n2;
  float d1;
  float d2;
  float o1;
  float o2;
  float *vmod;
  float *dmod;
  // parameters for modeling
  float *incomp;
  float *buoy;
};
/*^*/

#endif

static float* build_extended_model_2d(float const *mod, long n1, long n2, int next)
{

  sf_warning(" Extend the 2d models with absorbing boundaries..");

  long n1ext = n1+2*next;
  long n2ext = n2+2*next;
  long nelem = n1ext*n2ext;
  float* extd = sf_floatalloc(nelem);

  // core
  for (long i2=next,j2=0; i2<next+n2; i2++,j2++){
    for (long i1=next,j1=0; i1<next+n1; i1++,j1++){
      float const v = mod[j1+j2*n1];
      extd[i1+i2*n1ext] = v;
    }
  }

  // extend laterally
  for (long i2=0;i2<next;i2++){
    for (long i1=0,j1=next; i1<n1; i1++,j1++){
      extd[j1+i2*n1ext] = extd[j1+next*n1ext];
      extd[j1+(n2ext-1-i2)*n1ext] = extd[j1+(n2ext-1-next)*n1ext];
    }
  }

  // extend up and down
  for (long i2=0;i2<n2ext;i2++){
    for (long i1=0; i1<next; i1++){
      extd[i1+i2*n1ext] = extd[next+i2*n1ext];
      extd[(n1ext-1-i1)+i2*n1ext] = extd[(n1ext-1-next)+i2*n1ext];
    }
  }

  return extd;

}

void prepare_model_2d(mod_struct_t* mod,
                      in_para_struct_t para,
                      sf_axis axvel[2], sf_axis axden[2],
                      sf_file Fvmod, sf_file Fdmod)
/*< Sets up the specified model cube >*/
{

  if ((sf_n(axvel[0])!=sf_n(axden[0])) || (sf_n(axvel[1])!=sf_n(axden[1])))
    sf_error("Inconsistent model dimensions!");

  mod->n1 = sf_n(axvel[0]);
  mod->n2 = sf_n(axvel[1]);
  mod->d1 = sf_d(axvel[0]);
  mod->d2 = sf_d(axvel[1]);
  mod->o1 = sf_o(axvel[0]);
  mod->o2 = sf_o(axvel[1]);

  long nelem = mod->n1*mod->n2;

  mod->vmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->vmod,nelem,Fvmod);

  double vave = 0;
  for (int i=0; i<nelem; i++)
    vave += mod->vmod[i];
  vave /= (nelem);
  sf_warning("Velocity Model average value = %g",vave);

  mod->dmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->dmod,nelem,Fdmod);

  double dave = 0;
  for (int i=0; i<nelem; i++)
    dave += mod->dmod[i];
  dave /= (nelem);
  sf_warning("Velocity Model average value = %g",dave);

  // modeling parameters
  long n1 = mod->n1;
  long n2 = mod->n2;
  int nabc = para.nb;
  mod->incomp = build_extended_model_2d(mod->vmod,n1,n2,nabc);
  mod->buoy = build_extended_model_2d(mod->dmod,n1,n2,nabc);

  // compute incompressibility and buoyancy
  long nelemext = (n1+2*nabc)*(n2+2*nabc);
  for (long i=0; i<nelemext; i++){
    float v = mod->incomp[i];
    float d = mod->buoy[i];
    mod->incomp[i] = v*v*d;
    mod->buoy[i] = 1./d;
  }

}

void clear_model_2d(mod_struct_t* mod)
/*< Free the model parameter cubes >*/
{
  free(mod->vmod);
  free(mod->dmod);
  free(mod->incomp);
  free(mod->buoy);
}

void prepare_acquisition_2d( acq_struct_t* acq,
                          sf_axis axsou[2], sf_axis axrec[2], sf_axis axwav[2],
                          sf_file Fsou, sf_file Frec, sf_file Fwav)
/*< Read the acquisition geometry from files >*/
{
  //
  if (sf_n(axsou[0])!=2)
    sf_error("Wrong number of coordinates in the source file!");

  // shot coordinates
  acq->ns = sf_n(axsou[1]);
  acq->scoord = sf_floatalloc(2*acq->ns);
  for (long isou=0; isou<acq->ns; isou++)
    sf_floatread(acq->scoord+2*isou,2,Fsou);

  //
  if (sf_n(axrec[0])!=2)
    sf_error("Wrong number of coordinates in the receiver file!");

  // receiver coordinates
  acq->nr = sf_n(axrec[1]);
  acq->rcoord = sf_floatalloc(2*acq->nr);
  for (long irec=0; irec<acq->nr; irec++)
    sf_floatread(acq->rcoord+2*irec,2,Frec);

  // wavelet parameters
  acq->nt = sf_n(axwav[1]);
  acq->dt = sf_d(axwav[1]);
  acq->ot = sf_d(axwav[1]);

  long nsouwav = sf_n(axwav[0]);
  long nwavsamp = nsouwav*acq->nt;

  if (nsouwav==1){
    sf_warning("Using the same wavelet for all shots!");
  }
  else{
    if (nsouwav==acq->ns){
      sf_warning("Every shot has a different wavelet");
    }
    else{
      sf_error("Inconsistent number of wavelets and shots");
      return;
    }
  }

  acq->wav = sf_floatalloc(nwavsamp);
  sf_floatread(acq->wav,nwavsamp,Fwav);

}

/*
 * Evaluates the Bessel function (from Dave Hale's KaiserWindow class in the Java Mines TK)
 */
static double ino(double x)
{
  double s = 1.0;
  double ds = 1.0;
  double d = 0.0;
  do {
    d += 2.0;
    ds *= (x*x)/(d*d);
    s += ds;
  } while (ds>s*DBL_EPSILON);
  return s;
}

/*
 * Sinc function
 */
static double sinc(double x){

  if (x==0.0)
    return 1;
  else
    return sin(M_PI*x)/(M_PI*x);
}

/*
 * Kaiser window with fixed scale parameter and radius
 */
static double kwin(double x, double xmax){
  double b = 6.31;
  double xx = x*x;
  double xxmax = xmax*xmax;
  double scale = 1./ino(b);

  return (xx<xxmax)? scale*ino(b*sqrt(1.- xx/xxmax)):0.;

}

void set_sr_interpolation_coeffs(acq_struct_t * const acq, wfl_struct_t const * wfl)
/*< interpolation coefficients for source injection and receiver extraction >*/
{
  long nsous = acq->ns;
  long nrecs = acq->nr;

  float o1 = wfl->simO1;
  float o2 = wfl->simO2;
  float d1 = wfl->d1;
  float d2 = wfl->d2;

  acq->hicksSou1 = sf_floatalloc(8*nsous);
  acq->hicksSou2 = sf_floatalloc(8*nsous);
  acq->hicksRcv1 = sf_floatalloc(8*nrecs);
  acq->hicksRcv2 = sf_floatalloc(8*nrecs);

  for (long isou=0; isou<nsous; isou++){
    float x1s = acq->scoord[2*isou+1]; // z coordinate
    float x2s = acq->scoord[2*isou  ]; // x coordinate
    long ix1s = (x1s-o1)/d1;
    long ix2s = (x2s-o2)/d2;
    float rem1 = (x1s - (ix1s*d1+o1))/d1;
    float rem2 = (x2s - (ix2s*d2+o2))/d2;
    for (int i=-3,ii=0; i<=4; i++,ii++){
      acq->hicksSou1[ii+isou*8] = sinc(i-rem1)*kwin(i-rem1,4.5);
      acq->hicksSou2[ii+isou*8] = sinc(i-rem2)*kwin(i-rem2,4.5);
    }
  }

}

void clear_acq_2d(acq_struct_t *acq)
/*< Free the arrays in the acquisition structure >*/
{

  free(acq->hicksSou1);
  free(acq->hicksSou2);

  free(acq->scoord);
  free(acq->rcoord);

  free(acq->wav);
}

void prepare_wfl_2d(wfl_struct_t *wfl,mod_struct_t *mod, sf_file Fdata, sf_file Fwfl, in_para_struct_t para)
/*< Allocate the wavefield structure >*/
{

  wfl->modN1 = mod->n1;
  wfl->modN2 = mod->n2;

  wfl->nabc = para.nb;
  wfl->simN1 = mod->n1 + 2*para.nb;
  wfl->simN2 = mod->n2 + 2*para.nb;

  wfl->d1 = mod->d1;
  wfl->d2 = mod->d2;

  wfl->modO1 = mod->o1;
  wfl->modO2 = mod->o2;
  wfl->simO1 = mod->o1 - para.nb*mod->d1;
  wfl->simO2 = mod->o2 - para.nb*mod->d2;

  long nelem = wfl->simN1*wfl->simN2;
  long wflsize = nelem*sizeof(float);
  wfl->pc = sf_floatalloc(nelem);
  wfl->pp = sf_floatalloc(nelem);
  wfl->pa = sf_floatalloc(nelem);
  memset(wfl->pc,0,wflsize);
  memset(wfl->pp,0,wflsize);
  memset(wfl->pa,0,wflsize);

  wfl->v1c = sf_floatalloc(nelem);
  wfl->v1p = sf_floatalloc(nelem);
  wfl->v1a = sf_floatalloc(nelem);
  memset(wfl->v1c,0,wflsize);
  memset(wfl->v1p,0,wflsize);
  memset(wfl->v1a,0,wflsize);

  wfl->v2c = sf_floatalloc(nelem);
  wfl->v2p = sf_floatalloc(nelem);
  wfl->v2a = sf_floatalloc(nelem);
  memset(wfl->v2c,0,wflsize);
  memset(wfl->v2p,0,wflsize);
  memset(wfl->v2a,0,wflsize);

  // sponge
  wfl->tap1 = sf_floatalloc(wfl->simN1);
  wfl->tap2 = sf_floatalloc(wfl->simN2);

  wfl->bwfl = sf_floatalloc(wfl->modN1*wfl->modN2);

  wfl->Fdata = Fdata;
  wfl->Fwfl  = Fwfl;

}

void clear_wfl_2d(wfl_struct_t *wfl)
/*< Clear the wavefield structure >*/
{

  free(wfl->pc);
  free(wfl->pp);
  free(wfl->pa);

  free(wfl->v1c);
  free(wfl->v1p);
  free(wfl->v1a);

  free(wfl->v2c);
  free(wfl->v2p);
  free(wfl->v2a);

  free(wfl->tap1);
  free(wfl->tap2);

  free(wfl->bwfl);

}
